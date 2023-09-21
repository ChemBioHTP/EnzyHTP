"""TODO

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-08-17
"""
from typing import Union, List, Dict

from numpy.typing import ArrayLike
import numpy as np
import pandas as pd

from enzy_htp.core import _LOGGER, file_system as fs
from enzy_htp.structure import Structure, Atom

from enzy_htp import interface, config


def binding_energy(fname: str,
                   probe: str,
                   receptor: str,
                   probe_charge: int = None,
                   receptor_charge: int = None,
                   engine: str = 'xtb',
                   work_dir: str = None,
                   save_temp_files: bool = False,
                   **kwargs) -> float:
    """
    """

    if work_dir is None:
        work_dir = config['system.SCRATCH_DIR']

    #TODO(CJ): check here if the supplied method is supported

    return BINDING_ENERGY_ENGINE[engine](fname, probe, receptor, probe_charge, receptor_charge, work_dir, save_temp_files, **kwargs)


def xtb_binding_energy(fname: str, probe: str, receptor: str, probe_charge: int, receptor_charge: int, work_dir: str, save_temp_files: bool,
                       **kwargs) -> float:
    """

    """

    probe_xyz: str = f"{work_dir}/probe.xyz"
    receptor_xyz: str = f"{work_dir}/receptor.xyz"
    union_xyz: str = f"{work_dir}/union.xyz"
    probe_sdf: str = f"{work_dir}/probe.sdf"
    receptor_sdf: str = f"{work_dir}/receptor.sdf"
    union_xyz: str = f"{work_dir}/union.xyz"

    session = interface.pymol.new_session()

    interface.pymol.general_cmd(session, [('delete', 'all')])
    interface.pymol.create_cluster(session, fname, probe, outfile=probe_sdf, cap_strategy='CH3', work_dir=work_dir)

    if probe_charge is None:
        probe_charge: int = interface.bcl.calculate_formal_charge(probe_sdf)

    interface.pymol.general_cmd(session, [('delete', 'all')])
    interface.pymol.convert(session, probe_sdf, file_2=probe_xyz)

    interface.pymol.general_cmd(session, [('delete', 'all')])
    interface.pymol.create_cluster(session, fname, receptor, outfile=receptor_sdf, cap_strategy='CH3', work_dir=work_dir)

    if receptor_charge is None:
        receptor_charge: int = interface.bcl.calculate_formal_charge(receptor_sdf)

    interface.pymol.general_cmd(session, [('delete', 'all')])
    interface.pymol.convert(session, receptor_sdf, file_2=receptor_xyz)

    interface.pymol.general_cmd(session, [('delete', 'all')])
    interface.pymol.create_cluster(session, fname, f"({probe}) or ({receptor})", outfile=union_xyz, cap_strategy='CH3', work_dir=work_dir)

    union_charge = probe_charge + receptor_charge

    original_n_iter: int = config['xtb.N_ITER']

    union_energy: float = interface.xtb.single_point(union_xyz, charge=union_charge)
    probe_energy: float = interface.xtb.single_point(probe_xyz, charge=probe_charge)
    receptor_energy: float = interface.xtb.single_point(receptor_xyz, charge=receptor_charge)

    if not save_temp_files:
        fs.safe_rm(probe_xyz)
        fs.safe_rm(receptor_xyz)
        fs.safe_rm(union_xyz)
        fs.safe_rm(probe_sdf)
        fs.safe_rm(receptor_sdf)
        fs.safe_rm(union_xyz)

    return union_energy - (probe_energy + receptor_energy)


BINDING_ENERGY_ENGINE: Dict = {
    "xtb": xtb_binding_energy,
}
"""TODO(CJ)"""
