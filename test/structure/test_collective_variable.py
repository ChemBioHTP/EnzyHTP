"""Unit tests for enzy_htp.structure.collective_variable.

Tests cover:
    - CVTargets construction, iteration, serialization
    - DistanceCV / AngleCV / DihedralCV evaluation against known geometry
    - generate_window_targets (explicit list and range modes)
    - to_dict() / from_dict() round-trip for all CV types
    - AmberCV: to_structure_constraint, bind_topology, serialize_for_engine,
      _build_raw_rs_dict, _eval_amber_expr
    - PlumedCV: construction, serialize_for_engine, from_plumed, from_dict

Author: EnzyHTP
Date: 2026-04-29
"""
import os
import tempfile

import numpy as np
import pytest

from enzy_htp import PDBParser
from enzy_htp.structure.collective_variable import (
    CVTargets,
    CollectiveVariable,
    DistanceCV,
    AngleCV,
    DihedralCV,
    AmberCV,
    PlumedCV,
    _eval_amber_expr,
)

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = os.path.join(CURR_DIR, "data")

_PARSER = PDBParser()


# --------------------------------------------------------------------------
# Fixture: shared structure
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def ke_stru():
    """Load KE_07_R7_2_S.pdb once for the whole module."""
    pdb = os.path.join(DATA_DIR, "KE_07_R7_2_S.pdb")
    return _PARSER.get_structure(pdb)


# ==========================================================================
# CVTargets
# ==========================================================================

class TestCVTargets:
    def test_construction_from_list(self):
        t = CVTargets([1.0, 2.0, 3.0])
        assert list(t) == [1.0, 2.0, 3.0]
        assert len(t) == 3

    def test_getitem(self):
        t = CVTargets([10.0, 20.0])
        assert t[0] == 10.0
        assert t[-1] == 20.0

    def test_to_dict_from_dict_roundtrip(self):
        original = CVTargets([1.5, 2.5, 3.5])
        rebuilt = CVTargets.from_dict(original.to_dict())
        assert rebuilt.data == original.data

    def test_repr(self):
        t = CVTargets([1.0])
        assert "CVTargets" in repr(t)


# ==========================================================================
# generate_window_targets
# ==========================================================================

class TestGenerateWindowTargets:
    def test_explicit_list(self):
        cv = DistanceCV("B.254.H2", "A.101.OE2")
        tgts = cv.generate_window_targets([1.0, 2.0, 3.0])
        assert tgts.data == [1.0, 2.0, 3.0]

    def test_range_mode(self):
        cv = DistanceCV("B.254.H2", "A.101.OE2")
        tgts = cv.generate_window_targets(start=1.0, end=3.0, step=0.5)
        assert len(tgts) == 5
        assert np.isclose(tgts[0], 1.0)
        assert np.isclose(tgts[-1], 3.0)

    def test_range_mode_step_downward(self):
        cv = DistanceCV("B.254.H2", "A.101.OE2")
        tgts = cv.generate_window_targets(start=5.0, end=3.0, step=-1.0)
        assert len(tgts) == 3
        assert np.isclose(tgts[0], 5.0)
        assert np.isclose(tgts[-1], 3.0)

    def test_missing_args_raises(self):
        cv = DistanceCV("B.254.H2", "A.101.OE2")
        with pytest.raises(ValueError):
            cv.generate_window_targets()  # no args provided

    def test_partial_range_args_raises(self):
        cv = DistanceCV("B.254.H2", "A.101.OE2")
        with pytest.raises(ValueError):
            cv.generate_window_targets(start=1.0)  # missing end and step


# ==========================================================================
# DistanceCV
# ==========================================================================

class TestDistanceCV:
    def test_evaluate_on_structure(self, ke_stru):
        """Verified answer matches existing DistanceConstraint test."""
        cv = DistanceCV("B.254.H2", "A.101.OE2", name="d_test")
        val = cv.evaluate_on_structure(ke_stru)
        assert np.isclose(val, 2.0239901185529554, atol=1e-4)

    def test_unit(self):
        cv = DistanceCV("A.1.CA", "A.2.CA")
        assert cv.unit == "angstrom"

    def test_name(self):
        cv = DistanceCV("A.1.CA", "A.2.CA", name="bond_dist")
        assert cv.name == "bond_dist"

    def test_evaluate_on_frame_delegates(self, ke_stru):
        cv = DistanceCV("B.254.H2", "A.101.OE2")
        assert cv.evaluate_on_frame(ke_stru) == cv.evaluate_on_structure(ke_stru)

    def test_to_dict_from_dict_roundtrip(self):
        cv = DistanceCV("A.100.CA", "B.200.CB", name="myDist")
        d = cv.to_dict()
        assert d["cv_type"] == "distance"
        cv2 = CollectiveVariable.from_dict(d)
        assert isinstance(cv2, DistanceCV)
        assert cv2.atom_1_key == cv.atom_1_key
        assert cv2.atom_2_key == cv.atom_2_key
        assert cv2.name == cv.name

    def test_serialize_for_engine_non_amber(self):
        cv = DistanceCV("A.1.CA", "A.2.CA")
        result = cv.serialize_for_engine("plumed", "/tmp", 2.5)
        assert result["cv_type"] == "distance"
        assert result["window_target"] == 2.5

    def test_serialize_for_engine_amber_raises(self):
        cv = DistanceCV("A.1.CA", "A.2.CA")
        with pytest.raises(NotImplementedError):
            cv.serialize_for_engine("amber", "/tmp", 2.5)

    def test_missing_atom_key_raises(self, ke_stru):
        cv = DistanceCV("Z.999.XYZ", "A.101.OE2")
        with pytest.raises((KeyError, AttributeError, TypeError)):
            cv.evaluate_on_structure(ke_stru)


# ==========================================================================
# AngleCV
# ==========================================================================

class TestAngleCV:
    def test_evaluate_on_structure(self, ke_stru):
        """Angle B.254.H2–B.254.CAE–A.101.OE2 (via central atom B.254.H2)."""
        cv = AngleCV("B.254.CAE", "B.254.H2", "A.101.OE2", name="a_test")
        val = cv.evaluate_on_structure(ke_stru)
        # Central atom is B.254.H2; angle should be non-trivial
        assert 0.0 < val < 180.0

    def test_unit(self):
        cv = AngleCV("A.1.CA", "A.2.CA", "A.3.CA")
        assert cv.unit == "degree"

    def test_to_dict_from_dict_roundtrip(self):
        cv = AngleCV("A.1.CA", "A.2.CA", "A.3.CA", name="myAngle")
        d = cv.to_dict()
        cv2 = CollectiveVariable.from_dict(d)
        assert isinstance(cv2, AngleCV)
        assert cv2.atom_1_key == cv.atom_1_key
        assert cv2.atom_3_key == cv.atom_3_key
        assert cv2.name == cv.name


# ==========================================================================
# DihedralCV
# ==========================================================================

class TestDihedralCV:
    def test_evaluate_on_structure(self, ke_stru):
        """Answer matches the existing DihedralConstraint test."""
        cv = DihedralCV("B.254.CAE", "B.254.H2", "A.101.OE2", "A.101.CA")
        val = cv.evaluate_on_structure(ke_stru)
        assert np.isclose(val, 4.740006673137136, atol=1e-4)

    def test_unit(self):
        cv = DihedralCV("A.1.CA", "A.2.CA", "A.3.CA", "A.4.CA")
        assert cv.unit == "degree"

    def test_to_dict_from_dict_roundtrip(self):
        cv = DihedralCV("A.1.CA", "A.2.CB", "A.3.N", "A.4.C", name="myDih")
        d = cv.to_dict()
        cv2 = CollectiveVariable.from_dict(d)
        assert isinstance(cv2, DihedralCV)
        assert cv2.atom_4_key == cv.atom_4_key
        assert cv2.name == cv.name


# ==========================================================================
# _eval_amber_expr helper
# ==========================================================================

class TestEvalAmberExpr:
    def test_float_passthrough(self):
        assert _eval_amber_expr(3.0, 5.0) == 3.0

    def test_int_passthrough(self):
        assert _eval_amber_expr(4, 7.0) == 4.0

    def test_x_minus_offset(self):
        assert np.isclose(_eval_amber_expr("x-1", 5.0), 4.0)

    def test_x_plus_offset(self):
        assert np.isclose(_eval_amber_expr("x+0.5", 5.0), 5.5)

    def test_x_identity(self):
        assert np.isclose(_eval_amber_expr("x", 2.3), 2.3)

    def test_plain_string_float(self):
        assert _eval_amber_expr("3.14", 99.0) == pytest.approx(3.14)


# ==========================================================================
# AmberCV
# ==========================================================================

class TestAmberCV:
    @pytest.fixture
    def simple_amber_cv(self):
        dist_cv = DistanceCV("B.254.H2", "A.101.OE2", name="rxn_dist")
        return AmberCV(
            cv=dist_cv,
            amber_params={
                "rk2": 100.0,
                "rk3": 100.0,
                "r1": "x-1",
                "r2": "x",
                "r3": "x",
                "r4": "x+1",
            },
            name="amber_test",
        )

    def test_name(self, simple_amber_cv):
        assert simple_amber_cv.name == "amber_test"

    def test_unit_delegates_to_wrapped_cv(self, simple_amber_cv):
        assert simple_amber_cv.unit == "angstrom"

    def test_evaluate_on_structure_delegates(self, ke_stru, simple_amber_cv):
        val = simple_amber_cv.evaluate_on_structure(ke_stru)
        assert np.isclose(val, 2.0239901185529554, atol=1e-4)

    def test_to_dict_from_dict_roundtrip(self, simple_amber_cv):
        d = simple_amber_cv.to_dict()
        assert d["cv_type"] == "amber"
        cv2 = CollectiveVariable.from_dict(d)
        assert isinstance(cv2, AmberCV)
        assert cv2.name == simple_amber_cv.name
        assert isinstance(cv2.cv, DistanceCV)
        assert cv2.cv.atom_1_key == simple_amber_cv.cv.atom_1_key

    def test_serialize_for_engine_amber_writes_rs_file(self, simple_amber_cv):
        with tempfile.TemporaryDirectory() as tmpdir:
            payload = simple_amber_cv.serialize_for_engine("amber", tmpdir, 2.5)
            assert "disang" in payload
            assert os.path.exists(payload["disang"])
            with open(payload["disang"]) as f:
                content = f.read()
            assert "&rst" in content
            assert "rk2=100.0000" in content
            # r2 should equal x=2.5
            assert "r2=2.5000" in content
            # r1 should equal x-1=1.5
            assert "r1=1.5000" in content

    def test_serialize_for_engine_wrong_engine_raises(self, simple_amber_cv):
        with pytest.raises(NotImplementedError):
            simple_amber_cv.serialize_for_engine("gromacs", "/tmp", 2.5)

    def test_bind_topology_returns_new_instance(self, ke_stru, simple_amber_cv):
        bound = simple_amber_cv.bind_topology(ke_stru)
        assert bound is not simple_amber_cv
        assert "_resolved_iat" in bound.amber_params_

    def test_bind_topology_resolved_iat_are_ints(self, ke_stru, simple_amber_cv):
        bound = simple_amber_cv.bind_topology(ke_stru)
        iat = bound.amber_params_["_resolved_iat"]
        assert len(iat) == 2
        assert all(isinstance(idx, int) for idx in iat)

    def test_bound_serialize_uses_integer_iat(self, ke_stru, simple_amber_cv):
        bound = simple_amber_cv.bind_topology(ke_stru)
        with tempfile.TemporaryDirectory() as tmpdir:
            payload = bound.serialize_for_engine("amber", tmpdir, 2.5)
            with open(payload["disang"]) as f:
                content = f.read()
            # Should contain the resolved integer atom indices, not key strings
            iat = bound.amber_params_["_resolved_iat"]
            assert str(iat[0]) in content
            assert str(iat[1]) in content

    def test_to_structure_constraint_type(self, ke_stru, simple_amber_cv):
        from enzy_htp.structure.structure_constraint import DistanceConstraint
        cons = simple_amber_cv.to_structure_constraint(ke_stru, window_target=2.5)
        assert isinstance(cons, DistanceConstraint)

    def test_to_structure_constraint_target_value(self, ke_stru, simple_amber_cv):
        cons = simple_amber_cv.to_structure_constraint(ke_stru, window_target=3.7)
        assert cons.target_value == 3.7

    def test_to_structure_constraint_amber_params_propagated(self, ke_stru, simple_amber_cv):
        cons = simple_amber_cv.to_structure_constraint(ke_stru, window_target=2.5)
        assert "rk2" in cons.params["amber"]
        assert cons.params["amber"]["rk2"] == 100.0
        # r2 should be the expression 'x', not a resolved float
        assert cons.params["amber"]["r2"] == "x"

    def test_to_structure_constraint_angle_cv(self, ke_stru):
        from enzy_htp.structure.structure_constraint import AngleConstraint
        angle_cv = AngleCV("B.254.CAE", "B.254.H2", "A.101.OE2", name="a")
        amber_cv = AmberCV(
            cv=angle_cv,
            amber_params={"rk2": 50.0, "rk3": 50.0, "r2": "x"},
        )
        cons = amber_cv.to_structure_constraint(ke_stru, window_target=120.0)
        assert isinstance(cons, AngleConstraint)
        assert cons.target_value == 120.0

    def test_to_structure_constraint_dihedral_cv(self, ke_stru):
        from enzy_htp.structure.structure_constraint import DihedralConstraint
        dih_cv = DihedralCV("B.254.CAE", "B.254.H2", "A.101.OE2", "A.101.CA")
        amber_cv = AmberCV(cv=dih_cv, amber_params={"rk2": 50.0, "rk3": 50.0})
        cons = amber_cv.to_structure_constraint(ke_stru, window_target=0.0)
        assert isinstance(cons, DihedralConstraint)

    def test_generate_window_targets_range(self, simple_amber_cv):
        tgts = simple_amber_cv.generate_window_targets(start=1.0, end=3.0, step=0.5)
        assert len(tgts) == 5

    def test_unsupported_engine_raises(self, simple_amber_cv):
        with pytest.raises(NotImplementedError):
            simple_amber_cv.serialize_for_engine("plumed", "/tmp", 2.5)


# ==========================================================================
# PlumedCV
# ==========================================================================

class TestPlumedCV:
    PLUMED_STR = "d: DISTANCE ATOMS=1,10\nPRINT ARG=d FILE=colvar.dat"

    def test_construction(self):
        cv = PlumedCV(self.PLUMED_STR, name="plumed_dist", unit="nanometer")
        assert cv.name == "plumed_dist"
        assert cv.unit == "nanometer"
        assert cv.plumed_string == self.PLUMED_STR

    def test_label_inference(self):
        cv = PlumedCV("d: DISTANCE ATOMS=1,10")
        assert cv.label == "d"

    def test_label_inference_no_colon(self):
        cv = PlumedCV("DISTANCE ATOMS=1,10")
        assert cv.label is None

    def test_from_plumed_classmethod(self):
        cv = PlumedCV.from_plumed(self.PLUMED_STR, name="test")
        assert isinstance(cv, PlumedCV)
        assert cv.plumed_string == self.PLUMED_STR

    def test_evaluate_on_structure_raises(self, ke_stru):
        cv = PlumedCV(self.PLUMED_STR)
        with pytest.raises(NotImplementedError):
            cv.evaluate_on_structure(ke_stru)

    def test_serialize_for_engine_plumed_writes_dat(self):
        cv = PlumedCV(self.PLUMED_STR, name="p")
        with tempfile.TemporaryDirectory() as tmpdir:
            payload = cv.serialize_for_engine("plumed", tmpdir, 1.5)
            assert "plumed_dat" in payload
            assert os.path.exists(payload["plumed_dat"])
            with open(payload["plumed_dat"]) as f:
                content = f.read()
            assert "window_target = 1.5" in content
            assert "DISTANCE ATOMS=1,10" in content

    def test_serialize_for_engine_wrong_engine_raises(self):
        cv = PlumedCV(self.PLUMED_STR)
        with pytest.raises(NotImplementedError):
            cv.serialize_for_engine("amber", "/tmp", 1.5)

    def test_to_dict_from_dict_roundtrip(self):
        cv = PlumedCV(self.PLUMED_STR, name="plumed_cv_test", unit="nm", label="d")
        d = cv.to_dict()
        assert d["cv_type"] == "plumed"
        cv2 = CollectiveVariable.from_dict(d)
        assert isinstance(cv2, PlumedCV)
        assert cv2.plumed_string == cv.plumed_string
        assert cv2.name == cv.name
        assert cv2.unit == cv.unit
        assert cv2.label == cv.label


# ==========================================================================
# CollectiveVariable.from_dict dispatch
# ==========================================================================

class TestFromDictDispatch:
    def test_dispatch_distance(self):
        d = {"cv_type": "distance", "atom_1_key": "A.1.CA", "atom_2_key": "A.2.CA"}
        cv = CollectiveVariable.from_dict(d)
        assert isinstance(cv, DistanceCV)

    def test_dispatch_angle(self):
        d = {"cv_type": "angle", "atom_1_key": "A.1.CA", "atom_2_key": "A.2.CA",
             "atom_3_key": "A.3.CA"}
        cv = CollectiveVariable.from_dict(d)
        assert isinstance(cv, AngleCV)

    def test_dispatch_dihedral(self):
        d = {"cv_type": "dihedral", "atom_1_key": "A.1.CA", "atom_2_key": "A.2.CA",
             "atom_3_key": "A.3.CA", "atom_4_key": "A.4.CA"}
        cv = CollectiveVariable.from_dict(d)
        assert isinstance(cv, DihedralCV)

    def test_dispatch_amber(self):
        wrapped = {"cv_type": "distance", "atom_1_key": "A.1.CA", "atom_2_key": "A.2.CA"}
        d = {"cv_type": "amber", "wrapped_cv": wrapped, "amber_params": {"rk2": 100.0}}
        cv = CollectiveVariable.from_dict(d)
        assert isinstance(cv, AmberCV)

    def test_dispatch_plumed(self):
        d = {"cv_type": "plumed", "plumed_string": "d: DISTANCE ATOMS=1,2"}
        cv = CollectiveVariable.from_dict(d)
        assert isinstance(cv, PlumedCV)

    def test_unknown_cv_type_raises(self):
        with pytest.raises(KeyError):
            CollectiveVariable.from_dict({"cv_type": "nonexistent_cv_type"})
