"""Testing enzy_htp.geometry.sampling.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-09-25
"""

import glob
import pytest
import os
import numpy as np

from enzy_htp.core.clusters.accre import Accre
import enzy_htp.core.file_system as fs
from enzy_htp.core.general import load_obj
from enzy_htp.structure import structure_constraint as stru_cons
from enzy_htp.geometry import md_simulation, equi_md_sampling, deployable_equi_md_sampling, md_energy_injection, deployable_md_energy_injection
from enzy_htp import interface
from enzy_htp import PDBParser
from enzy_htp.structure import Atom, Residue, Chain, Structure, StruSelection
from enzy_htp._interface.amber_interface import AmberRestartParser
from enzy_htp._interface.handle_types import MolDynParameter, MolDynResult
from enzy_htp.core.exception import InconsistentMDEngine

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data/"
WORK_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/work_dir/"
sp = PDBParser()
amber_interface = interface.amber


class PicklableFakeParam(MolDynParameter):
    def __init__(self, inpcrd_path="initial.rst", prmtop_path="top.prmtop", ncaa_chrgspin_mapper=None):
        self._inpcrd = inpcrd_path
        self._prmtop = prmtop_path
        self.ncaa_chrgspin_mapper = ncaa_chrgspin_mapper or {}

    @property
    def engine(self):
        return "Amber"

    @property
    def topology_file(self):
        return self._prmtop

    @property
    def file_list(self):
        return [self._inpcrd, self._prmtop]

    @property
    def input_coordinate_file(self):
        return self._inpcrd


class PicklableFakeParamMethod:
    def __init__(self, temp_dir):
        self._parameterizer_temp_dir = str(temp_dir)

    @property
    def engine(self):
        return "Amber"

    @property
    def parameterizer_temp_dir(self):
        return self._parameterizer_temp_dir

    @parameterizer_temp_dir.setter
    def parameterizer_temp_dir(self, value):
        self._parameterizer_temp_dir = value

    @property
    def parent_interface(self):
        return amber_interface

    def run(self, stru):
        return PicklableFakeParam()


def _build_tiny_structure():
    atom_1 = Atom(name="C1", coord=(0.0, 0.0, 0.0), idx=1, element="C")
    atom_2 = Atom(name="O1", coord=(1.0, 0.0, 0.0), idx=2, element="O")
    atom_3 = Atom(name="H1", coord=(0.0, 1.0, 0.0), idx=3, element="H")
    residue = Residue(1, "LIG", [atom_1, atom_2, atom_3])
    chain = Chain("A", [residue])
    return Structure([chain])


def _patch_amber_atom_index(monkeypatch, mapper=None):
    if mapper is None:
        mapper = {}

    def get_amber_atom_index(atoms):
        return [mapper.get(atom.idx, atom.idx) for atom in atoms]

    monkeypatch.setattr(amber_interface, "get_amber_atom_index", get_amber_atom_index)


def _write_tiny_restart(path, title, coordinates, velocities=None):
    AmberRestartParser().write_restart(
        {
            "title": title,
            "atom_count": len(coordinates),
            "coordinates": coordinates,
            "velocities": velocities,
            "box": None,
            "time": None,
        },
        path,
    )


def _driven_region_kinetic_energy(selection, restart_path):
    restart = AmberRestartParser().get_restart(restart_path)
    assert restart["velocities"] is not None
    energy = 0.0
    amber_atom_indexes = amber_interface.get_amber_atom_index(selection.atoms)
    for atom, amber_atom_index in zip(selection.atoms, amber_atom_indexes):
        velocity = restart["velocities"][amber_atom_index - 1]
        energy += amber_interface._get_atom_mass(atom) * float(np.dot(velocity, velocity))
    return energy


def _assert_driven_region_kinetic_energy_increased(selection, segment_metadata, min_ratio=1.10):
    before = _driven_region_kinetic_energy(selection, segment_metadata["restart_in"])
    after = _driven_region_kinetic_energy(selection, segment_metadata["restart_out"])
    assert after > before * min_ratio


def test_amber_interface_velocity_scaling_to_restart(monkeypatch):
    test_stru = _build_tiny_structure()
    selection = StruSelection(test_stru.chains[0].residues[0].atoms[:2])
    _patch_amber_atom_index(monkeypatch)
    restart_in = f"{WORK_DIR}/scale_input.rst"
    restart_out = f"{WORK_DIR}/scale_output.rst"
    initial_velocities = [
        [0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.1, 0.0, 0.0],
    ]

    _write_tiny_restart(
        restart_in,
        "scale",
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        initial_velocities,
    )
    metadata = amber_interface.perturb_restart_velocities(
        selection,
        restart_in,
        restart_out,
        target_temperature=1200.0,
        mode="velocity_scale",
        remove_drift=False,
    )
    scaled_restart = AmberRestartParser().get_restart(restart_out)

    assert metadata["scale_factor"] > 0
    assert np.allclose(scaled_restart["coordinates"], [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    assert scaled_restart["velocities"][0][0] == pytest.approx(initial_velocities[0][0] * metadata["scale_factor"])
    assert scaled_restart["velocities"][1][0] == pytest.approx(initial_velocities[1][0] * metadata["scale_factor"])
    assert scaled_restart["velocities"][2][0] == pytest.approx(initial_velocities[2][0])

    fs.clean_temp_file_n_dir([restart_in, restart_out])


def test_amber_interface_velocity_scaling_single_atom_preserves_nonzero_velocity(monkeypatch):
    test_stru = _build_tiny_structure()
    selection = StruSelection([test_stru.chains[0].residues[0].atoms[0]])
    _patch_amber_atom_index(monkeypatch)
    restart_in = f"{WORK_DIR}/scale_single_atom_input.rst"
    restart_out = f"{WORK_DIR}/scale_single_atom_output.rst"

    _write_tiny_restart(
        restart_in,
        "scale_single",
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        [[0.5, 0.0, 0.0], [0.1, 0.0, 0.0], [0.1, 0.0, 0.0]],
    )

    amber_interface.perturb_restart_velocities(
        selection,
        restart_in,
        restart_out,
        target_temperature=1200.0,
        mode="velocity_scale",
    )
    scaled_restart = AmberRestartParser().get_restart(restart_out)

    assert abs(scaled_restart["velocities"][0][0]) > 0.0

    fs.clean_temp_file_n_dir([restart_in, restart_out])


def test_amber_interface_velocity_scaling_uses_amber_restart_order(monkeypatch):
    test_stru = _build_tiny_structure()
    selection = StruSelection(test_stru.chains[0].residues[0].atoms[:2])
    _patch_amber_atom_index(monkeypatch, {1: 3, 2: 1})
    restart_in = f"{WORK_DIR}/scale_mapped_input.rst"
    restart_out = f"{WORK_DIR}/scale_mapped_output.rst"
    initial_velocities = [
        [0.5, 0.0, 0.0],
        [0.2, 0.0, 0.0],
        [0.8, 0.0, 0.0],
    ]

    _write_tiny_restart(
        restart_in,
        "scale_mapped",
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        initial_velocities,
    )

    metadata = amber_interface.perturb_restart_velocities(
        selection,
        restart_in,
        restart_out,
        target_temperature=1200.0,
        mode="velocity_scale",
        remove_drift=False,
    )
    scaled_restart = AmberRestartParser().get_restart(restart_out)

    assert metadata["amber_atom_indexes"] == [3, 1]
    assert scaled_restart["velocities"][0][0] == pytest.approx(initial_velocities[0][0] * metadata["scale_factor"])
    assert scaled_restart["velocities"][1][0] == pytest.approx(initial_velocities[1][0])
    assert scaled_restart["velocities"][2][0] == pytest.approx(initial_velocities[2][0] * metadata["scale_factor"])

    fs.clean_temp_file_n_dir([restart_in, restart_out])


def test_amber_interface_maxwell_reassignment_to_restart(monkeypatch):
    test_stru = _build_tiny_structure()
    selection = StruSelection(test_stru.chains[0].residues[0].atoms[:2])
    _patch_amber_atom_index(monkeypatch)
    restart_in = f"{WORK_DIR}/maxwell_input.rst"
    restart_out = f"{WORK_DIR}/maxwell_output.rst"

    _write_tiny_restart(
        restart_in,
        "maxwell",
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        [[0.2, 0.0, 0.0], [0.1, 0.0, 0.0], [0.3, 0.0, 0.0]],
    )

    metadata = amber_interface.perturb_restart_velocities(
        selection,
        restart_in,
        restart_out,
        target_temperature=800.0,
        mode="maxwell_reassign",
    )
    reassigned_restart = AmberRestartParser().get_restart(restart_out)

    assert metadata["mode"] == "maxwell_reassign"
    assert np.allclose(reassigned_restart["coordinates"], [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    assert reassigned_restart["velocities"] is not None
    assert len(reassigned_restart["velocities"]) == 3
    assert reassigned_restart["velocities"][2][0] == pytest.approx(0.3)

    fs.clean_temp_file_n_dir([restart_in, restart_out])


def test_amber_interface_maxwell_single_atom_preserves_nonzero_velocity(monkeypatch):
    test_stru = _build_tiny_structure()
    selection = StruSelection([test_stru.chains[0].residues[0].atoms[0]])
    _patch_amber_atom_index(monkeypatch)
    restart_in = f"{WORK_DIR}/maxwell_single_atom_input.rst"
    restart_out = f"{WORK_DIR}/maxwell_single_atom_output.rst"

    _write_tiny_restart(
        restart_in,
        "maxwell_single",
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        [[0.2, 0.0, 0.0], [0.1, 0.0, 0.0], [0.3, 0.0, 0.0]],
    )

    amber_interface.perturb_restart_velocities(
        selection,
        restart_in,
        restart_out,
        target_temperature=800.0,
        mode="maxwell_reassign",
    )
    reassigned_restart = AmberRestartParser().get_restart(restart_out)

    assert abs(reassigned_restart["velocities"][0][0]) > 0.0

    fs.clean_temp_file_n_dir([restart_in, restart_out])


def test_amber_interface_maxwell_requires_velocities():
    test_stru = _build_tiny_structure()
    selection = StruSelection(test_stru.chains[0].residues[0].atoms[:2])
    restart_in = f"{WORK_DIR}/maxwell_input_no_vel.rst"
    restart_out = f"{WORK_DIR}/maxwell_output_no_vel.rst"

    _write_tiny_restart(
        restart_in,
        "maxwell",
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
    )

    with pytest.raises(ValueError):
        amber_interface.perturb_restart_velocities(
            selection,
            restart_in,
            restart_out,
            target_temperature=800.0,
            mode="maxwell_reassign",
        )

    fs.clean_temp_file_n_dir([restart_in, restart_out])


def test_amber_interface_maxwell_rejects_unsupported_element(monkeypatch):
    atom_1 = Atom(name="X1", coord=(0.0, 0.0, 0.0), idx=1, element="Xe")
    atom_2 = Atom(name="C1", coord=(1.0, 0.0, 0.0), idx=2, element="C")
    residue = Residue(1, "LIG", [atom_1, atom_2])
    chain = Chain("A", [residue])
    test_stru = Structure([chain])
    selection = StruSelection([atom_1])
    _patch_amber_atom_index(monkeypatch)
    restart_in = f"{WORK_DIR}/maxwell_input_bad_ele.rst"
    restart_out = f"{WORK_DIR}/maxwell_output_bad_ele.rst"

    _write_tiny_restart(
        restart_in,
        "maxwell",
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        [[0.5, 0.0, 0.0], [0.5, 0.0, 0.0]],
    )

    with pytest.raises(ValueError):
        amber_interface.perturb_restart_velocities(
            selection,
            restart_in,
            restart_out,
            target_temperature=800.0,
            mode="maxwell_reassign",
        )

    fs.clean_temp_file_n_dir([restart_in, restart_out])


def test_amber_interface_velocity_scaling_requires_velocities():
    test_stru = _build_tiny_structure()
    selection = StruSelection(test_stru.chains[0].residues[0].atoms[:2])
    restart_in = f"{WORK_DIR}/scale_input_no_vel.rst"
    restart_out = f"{WORK_DIR}/scale_output_no_vel.rst"

    _write_tiny_restart(
        restart_in,
        "scale",
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
    )

    with pytest.raises(ValueError):
        amber_interface.perturb_restart_velocities(
            selection,
            restart_in,
            restart_out,
            target_temperature=1200.0,
            mode="velocity_scale",
        )

    fs.clean_temp_file_n_dir([restart_in, restart_out])


def test_amber_interface_velocity_scaling_rejects_unsupported_element(monkeypatch):
    atom_1 = Atom(name="X1", coord=(0.0, 0.0, 0.0), idx=1, element="Xe")
    atom_2 = Atom(name="C1", coord=(1.0, 0.0, 0.0), idx=2, element="C")
    residue = Residue(1, "LIG", [atom_1, atom_2])
    chain = Chain("A", [residue])
    test_stru = Structure([chain])
    selection = StruSelection([atom_1])
    _patch_amber_atom_index(monkeypatch)
    restart_in = f"{WORK_DIR}/scale_input_bad_ele.rst"
    restart_out = f"{WORK_DIR}/scale_output_bad_ele.rst"

    _write_tiny_restart(
        restart_in,
        "scale",
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        [[0.5, 0.0, 0.0], [0.5, 0.0, 0.0]],
    )

    with pytest.raises(ValueError):
        amber_interface.perturb_restart_velocities(
            selection,
            restart_in,
            restart_out,
            target_temperature=1200.0,
            mode="velocity_scale",
        )

    fs.clean_temp_file_n_dir([restart_in, restart_out])


@pytest.mark.parametrize(
    "injection_mode,driven_region,driving_temperature",
    [
        ("velocity_scale", "resi 1", 1200.0),
        ("maxwell_reassign", "resi 1", 800.0),
    ],
)
def test_md_energy_injection_api_increases_driven_kinetic_energy(
        monkeypatch,
        tmp_path,
        injection_mode,
        driven_region,
        driving_temperature,
    ):
    test_stru = _build_tiny_structure()
    _patch_amber_atom_index(monkeypatch)

    coordinates = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    base_velocities = [[0.01, 0.0, 0.0], [0.02, 0.0, 0.0], [0.03, 0.0, 0.0]]

    class RestartResult:
        def __init__(self, last_frame_file):
            self.last_frame_file = last_frame_file
            self.source = "amber"

    class RestartWritingStep:
        def __init__(self, name, length):
            self.name = name
            self.length = length
            self.restart = False
            self.if_report = True
            self.record_period = 0.1
            self.ascii_rst = False
            self.work_dir = ""

        @property
        def engine(self):
            return "Amber"

        def run(self, input_data):
            source_restart = getattr(input_data, "input_coordinate_file", None)
            if self.name.startswith("prod"):
                velocities = base_velocities
            elif source_restart and os.path.exists(source_restart):
                restart = AmberRestartParser().get_restart(source_restart)
                velocities = restart["velocities"].tolist()
            else:
                velocities = base_velocities
            rst_path = os.path.join(self.work_dir, f"{self.name}.rst")
            _write_tiny_restart(rst_path, self.name, coordinates, velocities)
            return RestartResult(rst_path)

    fake_steps = (
        RestartWritingStep("min", 0.0),
        RestartWritingStep("heat", 0.1),
        RestartWritingStep("equi1", 0.1),
        RestartWritingStep("equi2", 0.1),
        RestartWritingStep("prod", 0.2),
    )

    def fake_process(**kwargs):
        return kwargs["param_method"], fake_steps

    monkeypatch.setattr("enzy_htp.geometry.sampling._process_equi_md_sampling_arguments", fake_process)
    monkeypatch.setattr("enzy_htp._interface.amber_interface.random.gauss", lambda mu, sigma: sigma)

    result = md_energy_injection(
        stru=test_stru,
        param_method=PicklableFakeParamMethod(tmp_path),
        work_dir=str(tmp_path / f"energy_injection_{injection_mode}"),
        prod_time=0.6,
        segment_time=0.2,
        driven_region=driven_region,
        injection_mode=injection_mode,
        driving_temperature=driving_temperature,
        stitched_trajectory=False,
        remove_drift=False,
    )

    replica = result["replicas"][0]
    assert len(replica["production"]) == 3
    assert len(replica["injections"]) == 2
    assert all(not atom.is_hydrogen() for atom in result["driven_region"].atoms)
    for segment_metadata in replica["injections"]:
        assert segment_metadata["mode"] == injection_mode
        assert segment_metadata["target_temperature"] == driving_temperature
        _assert_driven_region_kinetic_energy_increased(result["driven_region"], segment_metadata)


def test_md_energy_injection_argument_validation(tmp_path):
    test_stru = _build_tiny_structure()

    class FakeParamMethod:
        @property
        def engine(self):
            return "Amber"

        @property
        def parameterizer_temp_dir(self):
            return str(tmp_path)

        @parameterizer_temp_dir.setter
        def parameterizer_temp_dir(self, value):
            pass

        @property
        def parent_interface(self):
            return amber_interface

        def run(self, stru):
            raise Exception("should not run for validation failures")

    param_method = FakeParamMethod()

    with pytest.raises(ValueError):
        md_energy_injection(
            stru=test_stru,
            param_method=param_method,
            engine="gromacs",
            driven_region="resi 1",
            driving_temperature=600.0,
        )

    with pytest.raises(ValueError):
        md_energy_injection(
            stru=test_stru,
            param_method=param_method,
            injection_mode="not_a_mode",
            driven_region="resi 1",
            driving_temperature=600.0,
        )

    with pytest.raises(ValueError):
        md_energy_injection(
            stru=test_stru,
            param_method=param_method,
            driving_temperature=600.0,
        )

    with pytest.raises(ValueError):
        md_energy_injection(
            stru=test_stru,
            param_method=param_method,
            driven_region="resi 1",
        )

    class FakeNonAmberParamMethod(FakeParamMethod):
        @property
        def engine(self):
            return "OpenMM"

    with pytest.raises(InconsistentMDEngine):
        md_energy_injection(
            stru=test_stru,
            param_method=FakeNonAmberParamMethod(),
            driven_region="resi 1",
            driving_temperature=600.0,
        )


def test_md_energy_injection_local_driver(monkeypatch, tmp_path):
    test_stru = _build_tiny_structure()

    class FakeResult:
        def __init__(self, last_frame_file, traj_file, traj_log_file):
            self.last_frame_file = last_frame_file
            self.traj_file = traj_file
            self.traj_log_file = traj_log_file
            self.source = "amber"

    class FakeStep:
        def __init__(self, name, length):
            self.name = name
            self.length = length
            self.restart = False
            self.if_report = True
            self.record_period = 0.1
            self.ascii_rst = False
            self.work_dir = ""

        @property
        def engine(self):
            return "Amber"

        def run(self, input_data):
            return FakeResult(
                os.path.join(self.work_dir, f"{self.name}.rst"),
                os.path.join(self.work_dir, f"{self.name}.nc"),
                os.path.join(self.work_dir, f"{self.name}.out"),
            )

    fake_steps = (
        FakeStep("min", 0.0),
        FakeStep("heat", 0.1),
        FakeStep("equi1", 0.1),
        FakeStep("equi2", 0.1),
        FakeStep("prod", 0.2),
    )

    process_calls = []

    def fake_process(**kwargs):
        process_calls.append(kwargs)
        return kwargs["param_method"], fake_steps

    scaling_calls = []
    maxwell_calls = []
    combine_calls = []

    def fake_perturb(selection, restart_in, restart_out, target_temperature, mode="maxwell_reassign", remove_drift=True):
        if mode == "velocity_scale":
            scaling_calls.append((restart_in, restart_out, target_temperature, [atom.idx for atom in selection.atoms]))
            return {
                "restart_in": restart_in,
                "restart_out": restart_out,
                "current_temperature": 300.0,
                "target_temperature": target_temperature,
                "scale_factor": 2.0,
                "mode": "velocity_scale",
            }
        maxwell_calls.append((restart_in, restart_out, target_temperature, remove_drift))
        return {
            "restart_in": restart_in,
            "restart_out": restart_out,
            "target_temperature": target_temperature,
            "mode": "maxwell_reassign",
        }

    def fake_combine(traj_paths, topology_path, out_path, frame1_pdb_path=None, autoimage=True):
        combine_calls.append((traj_paths, topology_path, out_path, frame1_pdb_path, autoimage))
        fs.safe_mkdir(os.path.dirname(out_path))
        with open(out_path, "w") as of:
            of.write("combined")
        if frame1_pdb_path is not None:
            with open(frame1_pdb_path, "w") as of:
                of.write("frame1")

    monkeypatch.setattr("enzy_htp.geometry.sampling._process_equi_md_sampling_arguments", fake_process)
    monkeypatch.setattr(amber_interface, "perturb_restart_velocities", fake_perturb)
    monkeypatch.setattr(amber_interface, "combine_traj_segments", fake_combine)

    result = md_energy_injection(
        stru=test_stru,
        param_method=PicklableFakeParamMethod(tmp_path),
        work_dir=str(tmp_path / "energy_injection"),
        prod_time=0.4,
        segment_time=0.2,
        temperature=300.0,
        driven_region="resi 1",
        driving_temperature=600.0,
        injection_mode="velocity_scale",
        parallel_runs=1,
    )

    assert len(result["replicas"]) == 1
    assert len(result["replicas"][0]["equilibration"]) == 4
    assert len(result["replicas"][0]["production"]) == 2
    assert len(result["replicas"][0]["injections"]) == 1
    assert isinstance(result["replicas"][0]["equilibration"][0], MolDynResult)
    assert isinstance(result["replicas"][0]["production"][0], MolDynResult)
    assert scaling_calls[0][2] == 600.0
    assert scaling_calls[0][3] == [1, 2]
    assert not maxwell_calls
    assert process_calls[0]["prod_time"] == pytest.approx(0.2)
    assert process_calls[0]["record_period"] == pytest.approx(0.02)
    assert fake_steps[-1].ascii_rst is True
    assert fake_steps[-1].restart is True
    assert fake_steps[-1].record_period == pytest.approx(0.02)
    assert len(combine_calls) == 1
    assert combine_calls[0][4] is True
    assert len(combine_calls[0][0]) == 2
    assert combine_calls[0][3] is None
    assert os.path.exists(result["replicas"][0]["combined_traj_file"])
    assert "combined_frame1_pdb_file" not in result["replicas"][0]
    assert os.path.exists(result["replicas"][0]["checkpoint_path"])
    assert os.path.dirname(result["replicas"][0]["checkpoint_path"]) == str(tmp_path / "energy_injection")
    assert os.path.exists(os.path.join(str(tmp_path / "energy_injection"), "energy_injection_result.pickle"))

    saved_result = load_obj(os.path.join(str(tmp_path / "energy_injection"), "energy_injection_result.pickle"))
    assert saved_result["replicas"][0]["combined_traj_file"].endswith("prod_npt_combined.nc")
    assert isinstance(saved_result["replicas"][0]["equilibration"][0], dict)


def test_md_energy_injection_local_driver_maxwell_default(monkeypatch, tmp_path):
    test_stru = _build_tiny_structure()

    class FakeResult:
        def __init__(self, last_frame_file):
            self.last_frame_file = last_frame_file

    class FakeStep:
        def __init__(self, name, length):
            self.name = name
            self.length = length
            self.restart = False
            self.if_report = True
            self.record_period = 0.1
            self.work_dir = ""

        @property
        def engine(self):
            return "Amber"

        def run(self, input_data):
            return FakeResult(os.path.join(self.work_dir, f"{self.name}.rst"))

    fake_steps = (
        FakeStep("min", 0.0),
        FakeStep("heat", 0.1),
        FakeStep("equi1", 0.1),
        FakeStep("equi2", 0.1),
        FakeStep("prod", 0.2),
    )

    def fake_process(**kwargs):
        return kwargs["param_method"], fake_steps

    scaling_calls = []
    maxwell_calls = []

    def fake_perturb(selection, restart_in, restart_out, target_temperature, mode="maxwell_reassign", remove_drift=True):
        if mode == "velocity_scale":
            scaling_calls.append(((selection, restart_in, restart_out, target_temperature), {"remove_drift": remove_drift}))
            return {}
        maxwell_calls.append((restart_in, restart_out, target_temperature, remove_drift, [atom.idx for atom in selection.atoms]))
        return {
            "restart_in": restart_in,
            "restart_out": restart_out,
            "target_temperature": target_temperature,
            "mode": "maxwell_reassign",
        }

    monkeypatch.setattr("enzy_htp.geometry.sampling._process_equi_md_sampling_arguments", fake_process)
    monkeypatch.setattr(amber_interface, "perturb_restart_velocities", fake_perturb)

    result = md_energy_injection(
        stru=test_stru,
        param_method=PicklableFakeParamMethod(tmp_path),
        work_dir=str(tmp_path / "energy_injection_maxwell"),
        prod_time=0.4,
        segment_time=0.2,
        temperature=300.0,
        driven_region="resi 1",
        driving_temperature=900.0,
        remove_drift=False,
        parallel_runs=1,
        stitched_trajectory=False,
    )

    assert len(maxwell_calls) == 1
    assert not scaling_calls
    assert maxwell_calls[0][2] == 900.0
    assert maxwell_calls[0][3] is False
    assert maxwell_calls[0][4] == [1, 2]


def test_md_energy_injection_checkpoint_resume(monkeypatch, tmp_path):
    test_stru = _build_tiny_structure()
    call_log = []
    phase = {"fail_once": True}

    class FakeResult:
        def __init__(self, last_frame_file):
            self.last_frame_file = last_frame_file

    class FakeStep:
        def __init__(self, name, length):
            self.name = name
            self.length = length
            self.restart = False
            self.if_report = True
            self.record_period = 0.1
            self.work_dir = ""

        @property
        def engine(self):
            return "Amber"

        def run(self, input_data):
            call_log.append(self.name)
            if self.name.endswith("seg_000001") and phase["fail_once"]:
                phase["fail_once"] = False
                raise RuntimeError("simulated interruption")
            return FakeResult(os.path.join(self.work_dir, f"{self.name}.rst"))

    fake_steps = (
        FakeStep("min", 0.0),
        FakeStep("heat", 0.1),
        FakeStep("equi1", 0.1),
        FakeStep("equi2", 0.1),
        FakeStep("prod", 0.2),
    )

    def fake_process(**kwargs):
        return kwargs["param_method"], fake_steps

    def fake_perturb(selection, restart_in, restart_out, target_temperature, mode="maxwell_reassign", remove_drift=True):
        return {
            "restart_in": restart_in,
            "restart_out": restart_out,
            "current_temperature": 300.0,
            "target_temperature": target_temperature,
            "scale_factor": 2.0,
            "mode": "velocity_scale",
        }

    monkeypatch.setattr("enzy_htp.geometry.sampling._process_equi_md_sampling_arguments", fake_process)
    monkeypatch.setattr(amber_interface, "perturb_restart_velocities", fake_perturb)

    with pytest.raises(RuntimeError):
        md_energy_injection(
            stru=test_stru,
            param_method=PicklableFakeParamMethod(tmp_path),
            work_dir=str(tmp_path / "resume"),
            prod_time=0.4,
            segment_time=0.2,
            temperature=300.0,
            driven_region="resi 1",
            driving_temperature=600.0,
            injection_mode="velocity_scale",
            parallel_runs=1,
            stitched_trajectory=False,
        )

    first_phase_calls = list(call_log)
    result = md_energy_injection(
        stru=test_stru,
        param_method=PicklableFakeParamMethod(tmp_path),
        work_dir=str(tmp_path / "resume"),
        prod_time=0.4,
        segment_time=0.2,
        temperature=300.0,
        driven_region="resi 1",
        driving_temperature=600.0,
        injection_mode="velocity_scale",
        parallel_runs=1,
        stitched_trajectory=False,
    )

    resumed_calls = call_log[len(first_phase_calls):]
    assert first_phase_calls == ["min", "heat", "equi1", "equi2", "prod_seg_000000", "prod_seg_000001"]
    assert resumed_calls == ["prod_seg_000001"]
    assert len(result["replicas"][0]["production"]) == 2
    checkpoint = load_obj(result["replicas"][0]["checkpoint_path"])
    assert checkpoint["stage"] == "complete"
    assert checkpoint["next_segment_idx"] == 2


def test_md_energy_injection_local_driver_keeps_constraint_topology(monkeypatch, tmp_path):
    test_stru = _build_tiny_structure()
    test_constraint = stru_cons.create_distance_constraint(
        test_stru.atoms[0],
        test_stru.atoms[1],
        1.0,
    )


    class FakeResult:
        def __init__(self, last_frame_file):
            self.last_frame_file = last_frame_file

    class FakeStep:
        def __init__(self, name, length, constrain):
            self.name = name
            self.length = length
            self.restart = False
            self.if_report = True
            self.record_period = 0.1
            self.work_dir = ""
            self.constrain = constrain

        @property
        def engine(self):
            return "Amber"

        def run(self, input_data):
            for cons in self.constrain:
                assert cons.atoms[0].root() is test_stru
            return FakeResult(os.path.join(self.work_dir, f"{self.name}.rst"))

    fake_steps = (
        FakeStep("min", 0.0, [test_constraint]),
        FakeStep("heat", 0.1, [test_constraint]),
        FakeStep("equi1", 0.1, [test_constraint]),
        FakeStep("equi2", 0.1, [test_constraint]),
        FakeStep("prod", 0.2, [test_constraint]),
    )

    def fake_process(**kwargs):
        return kwargs["param_method"], fake_steps

    def fake_perturb(selection, restart_in, restart_out, target_temperature, mode="maxwell_reassign", remove_drift=True):
        return {
            "restart_in": restart_in,
            "restart_out": restart_out,
            "target_temperature": target_temperature,
            "mode": "maxwell_reassign",
        }

    monkeypatch.setattr("enzy_htp.geometry.sampling._process_equi_md_sampling_arguments", fake_process)
    monkeypatch.setattr(amber_interface, "perturb_restart_velocities", fake_perturb)

    result = md_energy_injection(
        stru=test_stru,
        param_method=PicklableFakeParamMethod(tmp_path),
        work_dir=str(tmp_path / "energy_injection_constraints"),
        prod_time=0.4,
        segment_time=0.2,
        temperature=300.0,
        driven_region="resi 1",
        driving_temperature=600.0,
        parallel_runs=1,
        prod_constrain=[test_constraint],
        stitched_trajectory=False,
    )

    assert len(result["replicas"][0]["production"]) == 2

def test_deployable_md_energy_injection():
    test_stru = _build_tiny_structure()

    cluster_job_config = {
        "cluster": Accre(),
        "res_keywords": {
            "core_type": "gpu",
            "account": "csb_gpu_acc",
            "partition": "turing",
            "node_cores": "1",
            "mem_per_core": "8G",
            "walltime": "3-00:00:00",
        },
    }

    result = deployable_md_energy_injection(
        stru=test_stru,
        param_method=PicklableFakeParamMethod(WORK_DIR),
        work_dir=f"{WORK_DIR}deployable_energy_injection",
        cluster_job_config=cluster_job_config,
        prod_time=0.4,
        segment_time=0.2,
        temperature=300.0,
        driven_region="resi 1",
        driving_temperature=600.0,
        parallel_runs=2,
    )

    assert len(result["job_list"]) == 2
    assert len(result["main_files"]) == 2
    assert len(result["kwargs_files"]) == 2
    assert len(result["result_files"]) == 2
    assert all(os.path.exists(main_file) for main_file in result["main_files"])
    assert all(os.path.exists(kwargs_file) for kwargs_file in result["kwargs_files"])
    assert result["result_files"][0].endswith("energy_injection_result.pickle")
    assert "python -u" in result["job_list"][0].sub_script_str
    assert "energy_injection_main.py" in result["job_list"][0].sub_script_str
    assert "energy_injection_replica_child_main.py" not in result["job_list"][0].sub_script_str
    assert "#SBATCH --gres=gpu:1" in result["job_list"][0].sub_script_str
    saved_kwargs = load_obj(result["kwargs_files"][0])
    assert "parallel_method" not in saved_kwargs["md_energy_injection_kwargs"]
    assert saved_kwargs["md_energy_injection_kwargs"]["parallel_runs"] == 1
    assert saved_kwargs["md_energy_injection_kwargs"]["driving_temperature"] == 600.0
    assert saved_kwargs["md_energy_injection_kwargs"]["result_fname"] == "energy_injection_result.pickle"
    assert os.path.isabs(saved_kwargs["md_energy_injection_kwargs"]["work_dir"])
    assert os.path.basename(saved_kwargs["md_energy_injection_kwargs"]["work_dir"]) == "rep_000000"

    fs.safe_rmdir(f"{WORK_DIR}deployable_energy_injection")


def test_deployable_md_energy_injection_requires_cluster_job_config():
    test_stru = _build_tiny_structure()

    class FakeParamMethod:
        @property
        def engine(self):
            return "Amber"

        @property
        def parameterizer_temp_dir(self):
            return WORK_DIR

        @parameterizer_temp_dir.setter
        def parameterizer_temp_dir(self, value):
            pass

        @property
        def parent_interface(self):
            return amber_interface

        def run(self, stru):
            raise Exception("should not run during deployable assembly")

    with pytest.raises(ValueError):
        deployable_md_energy_injection(
            stru=test_stru,
            param_method=FakeParamMethod(),
            work_dir=f"{WORK_DIR}deployable_energy_injection_missing_cfg",
            cluster_job_config=None,
            driven_region="resi 1",
            driving_temperature=600.0,
        )


@pytest.mark.accre
@pytest.mark.long
def test_md_energy_injection_lv1():
    """Test a real segmented energy-injection workflow with Maxwell reassignment."""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_param_method = amber_interface.build_md_parameterizer(
        ncaa_param_lib_path=f"{WORK_DIR}/ncaa_lib",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {"account" : "csb_gpu_acc",
                         "partition" : "turing"}
    }

    result = md_energy_injection(
        stru=test_stru,
        param_method=test_param_method,
        work_dir=f"{WORK_DIR}MD_ENERGY_INJECTION/",
        prod_time=0.1,
        segment_time=0.05,
        temperature=300.0,
        driven_region="resi 254",
        injection_mode="maxwell_reassign",
        driving_temperature=600.0,
        parallel_runs=1,
        cluster_job_config=cluster_job_config,
        record_period=0.01,
    )

    replica = result["replicas"][0]
    assert len(replica["equilibration"]) == 4
    assert len(replica["production"]) == 2
    assert len(replica["injections"]) == 1
    assert os.path.exists(replica["production"][0].traj_file)
    assert os.path.exists(replica["production"][0].traj_log_file)
    assert os.path.exists(replica["production"][0].last_frame_file)
    assert os.path.getsize(replica["production"][0].traj_file) > 0
    assert os.path.getsize(replica["production"][0].traj_log_file) > 0
    assert os.path.getsize(replica["production"][0].last_frame_file) > 0
    assert replica["production"][0].last_frame_parser(replica["production"][0].last_frame_file)
    for segment_metadata in replica["injections"]:
        _assert_driven_region_kinetic_energy_increased(result["driven_region"], segment_metadata)

    # clean up
    # fs.safe_rmdir(f"{WORK_DIR}MD_ENERGY_INJECTION/")
    # fs.clean_temp_file_n_dir([
    # ] + glob.glob("slurm-*.out")
    # + glob.glob("scratch/amber_parameterizer/*")
    # + glob.glob(f"{WORK_DIR}/ncaa_lib/H5J*"))


@pytest.mark.accre
@pytest.mark.long
def test_md_simulation_amber_lv1():
    """Test running a non-replica MD.
    Using Amber & Accre as an example engine
    level 1:
    - no replica
    - no constrain
    took around 2 min. use the internal checking in translate() to
    make sure each step is successfully finished."""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_param_method = amber_interface.build_md_parameterizer(
        ncaa_param_lib_path=f"{WORK_DIR}/ncaa_lib",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {"account" : "csb_gpu_acc",
                         "partition" : "turing"}
    }
    step_1  = amber_interface.build_md_step(
        name="min",
        minimize=True,
        length=2000, # cycle
        cluster_job_config=cluster_job_config,
        core_type="gpu",)

    step_2 = amber_interface.build_md_step(
        name="equi_npt",
        length=0.001, # ns
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=300.0,)

    step_3 = amber_interface.build_md_step(
        name="prod_npt",
        length=0.05, # ns
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=300.0,
        restart=True,
        if_report=True,
        record_period=0.0005,)

    md_result = md_simulation(stru=test_stru,
                  param_method=test_param_method,
                  steps=[step_1, step_2, step_3],
                  parallel_runs=1,
                  job_check_period=10)

    # TODO may be also check MD traj is reasonable?

    # clean up
    fs.clean_temp_file_n_dir([
    ] + glob.glob("slurm-*.out")
    + glob.glob("scratch/amber_parameterizer/*")
    + glob.glob("MD/rep_0/*out")
    + glob.glob("MD/rep_0/*nc")
    + glob.glob("MD/rep_0/*rst")
    + glob.glob(f"{WORK_DIR}/ncaa_lib/H5J*"))

@pytest.mark.accre
@pytest.mark.long
def test_md_simulation_amber_lv2():
    """Test running a non-replica MD.
    Using Amber & Accre as an example engine
    level 1:
    - no replica
    - backbone constrain"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_param_method = amber_interface.build_md_parameterizer(
        ncaa_param_lib_path=f"{WORK_DIR}/ncaa_lib",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {"account" : "csb_gpu_acc",
                         "partition" : "turing"}
    }

    constrain = [stru_cons.create_backbone_freeze(test_stru)]

    step_1  = amber_interface.build_md_step(
        name="min",
        minimize=True,
        length=2000, # cycle
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        constrain=constrain,)

    step_2 = amber_interface.build_md_step(
        name="equi_npt",
        length=0.001, # ns
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=300.0,
        constrain=constrain,)

    step_3 = amber_interface.build_md_step(
        name="prod_npt",
        length=0.05, # ns
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=300.0,
        restart=True,
        if_report=True,
        record_period=0.0005,
        constrain=constrain,)

    md_result = md_simulation(
        stru=test_stru,
        param_method=test_param_method,
        steps=[step_1, step_2, step_3],
        parallel_runs=1,
        job_check_period=10)

    # clean up
    fs.clean_temp_file_n_dir([
    ] + glob.glob("slurm-*.out")
    + glob.glob("scratch/amber_parameterizer/*")
    + glob.glob("MD/rep_0/*out")
    + glob.glob("MD/rep_0/*nc")
    + glob.glob("MD/rep_0/*rst")
    + glob.glob(f"{WORK_DIR}/ncaa_lib/H5J*"))


@pytest.mark.accre
@pytest.mark.long
def test_md_simulation_amber_lv3():
    """Test running a non-replica MD.
    Using Amber & Accre as an example engine
    level 1:
    - no replica
    - backbone constrain
    - distance and angle constrain"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_param_method = amber_interface.build_md_parameterizer(
        ncaa_param_lib_path=f"{WORK_DIR}/ncaa_lib",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {"account" : "csb_gpu_acc",
                         "partition" : "turing"}
    }

    constrain = [stru_cons.create_backbone_freeze(test_stru),
                 stru_cons.create_distance_constraint(
                    "B.254.H2", "A.101.OE2", 2.4, test_stru),
                 stru_cons.create_angle_constraint(
                    "B.254.CAE", "B.254.H2", "A.101.OE2", 180.0, test_stru),]

    step_1  = amber_interface.build_md_step(
        name="min",
        minimize=True,
        length=2000, # cycle
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        constrain=constrain,)

    step_2 = amber_interface.build_md_step(
        name="equi_npt",
        length=0.001, # ns
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=300.0,
        constrain=constrain,)

    step_3 = amber_interface.build_md_step(
        name="prod_npt",
        length=0.05, # ns
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=300.0,
        restart=True,
        if_report=True,
        record_period=0.0005,
        constrain=constrain,)

    md_result = md_simulation(
        stru=test_stru,
        param_method=test_param_method,
        steps=[step_1, step_2, step_3],
        parallel_runs=1,
        job_check_period=10)

    assert len(glob.glob("MD/rep_0/*in")) == 0

    # clean up
    fs.clean_temp_file_n_dir([
    ] + glob.glob("slurm-*.out")
    + glob.glob("scratch/amber_parameterizer/*")
    + glob.glob("MD/rep_0/*out")
    + glob.glob("MD/rep_0/*nc")
    + glob.glob("MD/rep_0/*rst")
    + glob.glob("MD/rep_0/*rs")
    + glob.glob(f"{WORK_DIR}/ncaa_lib/H5J*"))


@pytest.mark.accre
@pytest.mark.long
def test_md_simulation_amber_3_repeat():
    """Test running a 3-replica MD.
    Using Amber & Accre as an example engine"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_param_method = amber_interface.build_md_parameterizer(
        ncaa_param_lib_path=f"{WORK_DIR}/ncaa_lib",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {"account" : "csb_gpu_acc",
                         "partition" : "turing"}
    }

    constrain = [stru_cons.create_backbone_freeze(test_stru)]

    step_1  = amber_interface.build_md_step(
        name="min",
        minimize=True,
        length=2000, # cycle
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        constrain=constrain,)

    step_2 = amber_interface.build_md_step(
        name="equi_npt",
        length=0.001, # ns
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=300.0,
        constrain=constrain,)

    step_3 = amber_interface.build_md_step(
        name="prod_npt",
        length=0.05, # ns
        cluster_job_config=cluster_job_config,
        core_type="gpu",
        temperature=300.0,
        restart=True,
        if_report=True,
        record_period=0.0005,
        constrain=constrain,)

    md_result = md_simulation(
        stru=test_stru,
        param_method=test_param_method,
        steps=[step_1, step_2, step_3],
        parallel_runs=3,
        job_check_period=10)

    # clean up
    fs.clean_temp_file_n_dir([
    ] + glob.glob("slurm-*.out")
    + glob.glob("scratch/amber_parameterizer/*")
    + glob.glob("MD/rep_0/*out")
    + glob.glob("MD/rep_0/*nc")
    + glob.glob("MD/rep_0/*rst")
    + glob.glob("MD/rep_1/*out")
    + glob.glob("MD/rep_1/*nc")
    + glob.glob("MD/rep_1/*rst")
    + glob.glob("MD/rep_2/*out")
    + glob.glob("MD/rep_2/*nc")
    + glob.glob("MD/rep_2/*rst")
    + glob.glob(f"{WORK_DIR}/ncaa_lib/H5J*"))

@pytest.mark.accre
@pytest.mark.long
def test_equi_md_sampling_lv1():
    """test for equi_md_sampling
    level 1: no constraint"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_param_method = amber_interface.build_md_parameterizer(
        ncaa_param_lib_path=f"{WORK_DIR}/ncaa_lib",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {"account" : "csb_gpu_acc",
                         "partition" : "turing"}
    }
    md_result = equi_md_sampling(
        stru = test_stru,
        param_method = test_param_method,
        cluster_job_config = cluster_job_config,
        job_check_period=10,
        # shorter sim for test
        prod_time=0.5,
        record_period=0.05,
        work_dir=f"{WORK_DIR}MD/")

    # clean up
    # fs.safe_rmdir(f"{WORK_DIR}MD/")
    # fs.clean_temp_file_n_dir([
    # ] + glob.glob("slurm-*.out")
    # + glob.glob("scratch/amber_parameterizer/*")
    # + glob.glob(f"{WORK_DIR}/ncaa_lib/H5J*"))


@pytest.mark.accre
def test_equi_md_sampling_lv2():
    """test for equi_md_sampling
    level 2:
    - geom constrain"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_param_method = amber_interface.build_md_parameterizer(
        ncaa_param_lib_path=f"{WORK_DIR}/ncaa_lib",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {"account" : "csb_gpu_acc",
                         "partition" : "turing"}
    }
    constrain = [stru_cons.create_distance_constraint(
                    "B.254.H2", "A.101.OE2", 2.4, test_stru),
                 stru_cons.create_angle_constraint(
                    "B.254.CAE", "B.254.H2", "A.101.OE2", 180.0, test_stru),]
    md_result = equi_md_sampling(
        stru = test_stru,
        param_method = test_param_method,
        cluster_job_config = cluster_job_config,
        job_check_period=10,
        prod_constrain=constrain,
        # shorter sim for test
        prod_time=0.5,
        record_period=0.05,
        work_dir=f"{WORK_DIR}MD/")

    # clean up
    fs.safe_rmdir(f"{WORK_DIR}MD/")
    fs.clean_temp_file_n_dir([
    ] + glob.glob("slurm-*.out")
    + glob.glob("scratch/amber_parameterizer/*")
    + glob.glob(f"{WORK_DIR}/ncaa_lib/H5J*"))

def test_equi_md_sampling_wrong_cons():
    """test for equi_md_sampling that uses wrong constraints"""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru_2 = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S_mut.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_param_method = amber_interface.build_md_parameterizer(
        ncaa_param_lib_path=f"{WORK_DIR}/ncaa_lib",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {"account" : "csb_gpu_acc",
                         "partition" : "turing"}
    }
    constrain = [stru_cons.create_distance_constraint(
                    "B.254.H2", "A.101.OE2", 2.4, test_stru_2),
                 stru_cons.create_angle_constraint(
                    "B.254.CAE", "B.254.H2", "A.101.OE2", 180.0, test_stru_2),]

    with pytest.raises(ValueError) as e:
        md_result = equi_md_sampling(
            stru = test_stru,
            param_method = test_param_method,
            cluster_job_config = cluster_job_config,
            job_check_period=10,
            prod_constrain=constrain,
            # shorter sim for test
            prod_time=0.5,
            record_period=0.05,
            work_dir=f"{WORK_DIR}MD/")

@pytest.mark.accre
def test_deployable_equi_md_sampling():
    """test for deployable_equi_md_sampling. simply make sure no error."""
    test_stru = sp.get_structure(f"{DATA_DIR}KE_07_R7_2_S.pdb")
    test_stru.assign_ncaa_chargespin({"H5J" : (0,1)})
    test_param_method = amber_interface.build_md_parameterizer(
        ncaa_param_lib_path=f"{WORK_DIR}/ncaa_lib",
    )
    cluster_job_config = {
        "cluster" : Accre(),
        "res_keywords" : {"account" : "csb_gpu_acc",
                         "partition" : "turing"}
    }
    constrain = [stru_cons.create_distance_constraint(
                    "B.254.H2", "A.101.OE2", 2.4, test_stru),
                 stru_cons.create_angle_constraint(
                    "B.254.CAE", "B.254.H2", "A.101.OE2", 180.0, test_stru),]
    md_result = deployable_equi_md_sampling(
        stru = test_stru,
        param_method = test_param_method,
        cluster_job_config = cluster_job_config,
        prod_constrain=constrain,
        # shorter sim for test
        prod_time=0.5,
        record_period=0.05,
        work_dir=f"{WORK_DIR}MD/")

    # clean up
    fs.safe_rmdir(f"{WORK_DIR}MD/")
    fs.clean_temp_file_n_dir(glob.glob(f"{WORK_DIR}/ncaa_lib/H5J*"))
