"""Defines the CollectiveVariable (CV) abstraction and concrete implementations for
EnzyHTP.  CVs are the reaction-coordinate definitions used in umbrella sampling and
related free-energy calculations.

Classes exported by this module:
    CVTargets           – ordered list of per-window target values
    CollectiveVariable  – abstract base class
    DistanceCV          – distance between two atom selections (Å)
    AngleCV             – angle between three atom selections (degrees)
    DihedralCV          – torsion angle between four atom selections (degrees)
    AmberCV             – CV with Amber DISANG / NMROPT restraint serialization
    PlumedCV            – CV declared via a raw PLUMED string

Author: EnzyHTP
Date: 2026-04-29
"""
from __future__ import annotations

import os
import math
from abc import ABC, abstractmethod
from copy import deepcopy
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from enzy_htp.core.logger import _LOGGER

# ---------------------------------------------------------------------------
# CVTargets
# ---------------------------------------------------------------------------

class CVTargets:
    """Ordered sequence of window target values for a CollectiveVariable.

    Attributes:
        data: list of float target values, one per umbrella window.
    """

    def __init__(self, data: List[float]) -> None:
        self.data: List[float] = list(data)

    # --- dunder helpers -------------------------------------------------
    def __len__(self) -> int:
        return len(self.data)

    def __iter__(self):
        return iter(self.data)

    def __getitem__(self, idx):
        return self.data[idx]

    def __repr__(self) -> str:
        return f"CVTargets(data={self.data!r})"

    # --- serialization --------------------------------------------------
    def to_dict(self) -> Dict[str, Any]:
        """Serialize to a plain dict."""
        return {"data": list(self.data)}

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "CVTargets":
        """Reconstruct from a plain dict produced by :py:meth:`to_dict`."""
        return cls(d["data"])


# ---------------------------------------------------------------------------
# CollectiveVariable (abstract base class)
# ---------------------------------------------------------------------------

# Registry maps cv_type string → subclass, populated by each subclass.
_CV_REGISTRY: Dict[str, type] = {}


def _register_cv(cls):
    """Class decorator that registers a CV subclass by its ``cv_type`` attribute."""
    _CV_REGISTRY[cls.cv_type] = cls
    return cls


class CollectiveVariable(ABC):
    """Abstract base class that all Collective Variable implementations inherit from.

    Subclasses must define:
        cv_type     – class-level str used for serialization dispatch
        name        – human-readable label (property or attribute)
        unit        – unit string, e.g. 'angstrom' or 'degree'
        evaluate_on_structure(structure) → float
        serialize_for_engine(engine, work_dir, window_target) → dict
        to_dict() → dict
        from_dict(data) → CollectiveVariable  (classmethod)

    The ``generate_window_targets`` helper is implemented in the base class and
    works for any concrete subclass.
    """

    cv_type: str = "base"  # overridden in subclasses

    # --- core properties (may be set as attributes or overridden as properties)

    @property
    @abstractmethod
    def name(self) -> str:
        """Human-readable label for this CV."""

    @property
    @abstractmethod
    def unit(self) -> str:
        """Unit string (e.g. 'angstrom', 'degree')."""

    # --- window target generation ---------------------------------------

    def generate_window_targets(
        self,
        data: Optional[Sequence[float]] = None,
        *,
        start: Optional[float] = None,
        end: Optional[float] = None,
        step: Optional[float] = None,
    ) -> CVTargets:
        """Generate a :class:`CVTargets` list for umbrella sampling windows.

        Two modes are supported:
            1. **Explicit list** – pass *data* as a sequence of floats.
            2. **Range** – pass *start*, *end*, and *step* (numpy-linspace
               style, i.e. both ends are inclusive). The number of points is
               ``round((end - start) / step) + 1``.

        Args:
            data:  Optional explicit list of target values.
            start: Start of range (inclusive).
            end:   End of range (inclusive).
            step:  Step between successive windows.

        Returns:
            A :class:`CVTargets` instance.

        Raises:
            ValueError: If neither *data* nor a complete (*start*, *end*, *step*)
                triple is provided.
        """
        if data is not None:
            return CVTargets(list(data))

        if start is not None and end is not None and step is not None:
            n = round(abs(end - start) / abs(step)) + 1
            targets = [start + i * step for i in range(n)]
            return CVTargets(targets)

        raise ValueError(
            "generate_window_targets() requires either `data` or all three of "
            "`start`, `end`, `step`."
        )

    # --- evaluation -------------------------------------------------

    @abstractmethod
    def evaluate_on_structure(self, structure) -> float:
        """Compute the CV scalar for the given :class:`~enzy_htp.structure.Structure`.

        Args:
            structure: A :class:`~enzy_htp.structure.Structure` object.

        Returns:
            The scalar value of this CV in ``self.unit``.
        """

    def evaluate_on_frame(self, frame) -> float:
        """Evaluate CV on a trajectory frame.

        By default delegates to :py:meth:`evaluate_on_structure`.  Subclasses
        may override to handle raw coordinate arrays.

        Args:
            frame: A :class:`~enzy_htp.structure.Structure` representing one
                trajectory frame.

        Returns:
            Scalar CV value.
        """
        return self.evaluate_on_structure(frame)

    # --- serialization --------------------------------------------------

    @abstractmethod
    def serialize_for_engine(
        self,
        engine: str,
        work_dir: str,
        window_target: float,
    ) -> Dict[str, Any]:
        """Produce engine-specific files/metadata for a given window target.

        Args:
            engine:         Target MD engine identifier, e.g. ``'amber'`` or
                            ``'plumed'``.
            work_dir:       Directory where engine-specific files should be
                            written.
            window_target:  The RC target value for this window (in ``self.unit``).

        Returns:
            A dict with engine-specific payload keys.
        """

    @abstractmethod
    def to_dict(self) -> Dict[str, Any]:
        """Serialize to a plain dict that fully describes this CV.

        The returned dict must include a ``cv_type`` key matching ``cls.cv_type``
        so that :py:meth:`from_dict` can reconstruct it.
        """

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "CollectiveVariable":
        """Reconstruct a :class:`CollectiveVariable` from a dict produced by
        :py:meth:`to_dict`.

        The ``cv_type`` key in *data* is used to dispatch to the correct
        subclass. This method may be called on the base class or on any
        concrete subclass.

        Args:
            data: Dict with at minimum a ``cv_type`` key.

        Returns:
            The reconstructed :class:`CollectiveVariable` instance.

        Raises:
            KeyError: If ``cv_type`` is unrecognised.
        """
        cv_type = data.get("cv_type", cls.cv_type)
        target_cls = _CV_REGISTRY.get(cv_type)
        if target_cls is None:
            raise KeyError(
                f"Unknown cv_type '{cv_type}'. "
                f"Available: {sorted(_CV_REGISTRY.keys())}"
            )
        return target_cls._from_dict_impl(data)

    @classmethod
    def _from_dict_impl(cls, data: Dict[str, Any]) -> "CollectiveVariable":
        """Subclass hook called by :py:meth:`from_dict`.  Override in each
        concrete subclass."""
        raise NotImplementedError(f"{cls.__name__}._from_dict_impl() not implemented.")

    # --- helpers --------------------------------------------------------

    def __repr__(self) -> str:
        return f"{type(self).__name__}(name={self.name!r}, unit={self.unit!r})"


# ---------------------------------------------------------------------------
# Helpers for atom-key–based CVs
# ---------------------------------------------------------------------------

def _resolve_atom(key: str, structure):
    """Return the :class:`~enzy_htp.structure.Atom` corresponding to *key* in
    *structure*.

    Args:
        key:       Atom key string understood by ``Structure.get()``,
                   e.g. ``'A.100.CA'``.
        structure: A :class:`~enzy_htp.structure.Structure` object.

    Returns:
        The matching :class:`~enzy_htp.structure.Atom`.

    Raises:
        KeyError: If *key* is not found in *structure*.
    """
    atom = structure.get(key)
    if atom is None:
        raise KeyError(f"Atom key '{key}' not found in structure.")
    return atom


# ---------------------------------------------------------------------------
# DistanceCV
# ---------------------------------------------------------------------------

@_register_cv
class DistanceCV(CollectiveVariable):
    """Distance between two atoms, measured in Ångströms.

    Atoms are referenced by their key strings (e.g. ``'A.100.CA'``) so the CV
    remains portable across different :class:`~enzy_htp.structure.Structure`
    instances that share the same topology.

    Attributes:
        name_:        Human-readable label.
        atom_1_key:   Key string for the first atom.
        atom_2_key:   Key string for the second atom.
    """

    cv_type: str = "distance"

    def __init__(
        self,
        atom_1_key: str,
        atom_2_key: str,
        name: str = "distance_cv",
    ) -> None:
        self.atom_1_key = atom_1_key
        self.atom_2_key = atom_2_key
        self.name_: str = name

    # --- CollectiveVariable interface -----------------------------------

    @property
    def name(self) -> str:
        return self.name_

    @property
    def unit(self) -> str:
        return "angstrom"

    def evaluate_on_structure(self, structure) -> float:
        """Compute inter-atomic distance in Å."""
        a1 = _resolve_atom(self.atom_1_key, structure)
        a2 = _resolve_atom(self.atom_2_key, structure)
        return float(a1.distance_to(a2))

    def serialize_for_engine(
        self,
        engine: str,
        work_dir: str,
        window_target: float,
    ) -> Dict[str, Any]:
        """Return metadata describing this CV for *engine*.

        For ``'amber'`` the caller should use :class:`AmberCV` which wraps this
        CV with the necessary force-constant parameters.  Calling
        ``serialize_for_engine`` directly on a ``DistanceCV`` with
        ``engine='amber'`` raises a ``NotImplementedError``; use
        :class:`AmberCV` instead.

        Args:
            engine:        Target engine string.
            work_dir:      Working directory (not used for distance metadata).
            window_target: Target value for this window.

        Returns:
            Dict with keys ``cv_type``, ``atom_1``, ``atom_2``,
            ``window_target``, and ``unit``.

        Raises:
            NotImplementedError: If *engine* is ``'amber'``.  Use
                :class:`AmberCV` instead.
        """
        if engine == "amber":
            raise NotImplementedError(
                "Direct Amber serialization of DistanceCV is not supported. "
                "Wrap this CV with AmberCV to emit DISANG restraint files."
            )
        return {
            "cv_type": self.cv_type,
            "atom_1": self.atom_1_key,
            "atom_2": self.atom_2_key,
            "window_target": window_target,
            "unit": self.unit,
        }

    def to_dict(self) -> Dict[str, Any]:
        return {
            "cv_type": self.cv_type,
            "name": self.name_,
            "atom_1_key": self.atom_1_key,
            "atom_2_key": self.atom_2_key,
        }

    @classmethod
    def _from_dict_impl(cls, data: Dict[str, Any]) -> "DistanceCV":
        return cls(
            atom_1_key=data["atom_1_key"],
            atom_2_key=data["atom_2_key"],
            name=data.get("name", "distance_cv"),
        )


# ---------------------------------------------------------------------------
# AngleCV
# ---------------------------------------------------------------------------

@_register_cv
class AngleCV(CollectiveVariable):
    """Angle formed by three atoms, measured in degrees.

    The angle is computed at the central atom (atom_2).

    Attributes:
        name_:        Human-readable label.
        atom_1_key:   Key string for the first atom.
        atom_2_key:   Key string for the central atom.
        atom_3_key:   Key string for the third atom.
    """

    cv_type: str = "angle"

    def __init__(
        self,
        atom_1_key: str,
        atom_2_key: str,
        atom_3_key: str,
        name: str = "angle_cv",
    ) -> None:
        self.atom_1_key = atom_1_key
        self.atom_2_key = atom_2_key
        self.atom_3_key = atom_3_key
        self.name_: str = name

    @property
    def name(self) -> str:
        return self.name_

    @property
    def unit(self) -> str:
        return "degree"

    def evaluate_on_structure(self, structure) -> float:
        """Compute the angle in degrees using the dot-product formula."""
        a1 = _resolve_atom(self.atom_1_key, structure)
        a2 = _resolve_atom(self.atom_2_key, structure)
        a3 = _resolve_atom(self.atom_3_key, structure)
        # Atom.angle_with(point_1, point_2) computes angle self–point_1–point_2
        # We want a1–a2–a3, so we use a2.angle_with(a1, a3)
        return float(a2.angle_with(a1, a3))

    def serialize_for_engine(
        self,
        engine: str,
        work_dir: str,
        window_target: float,
    ) -> Dict[str, Any]:
        """Return metadata for *engine*.  Direct Amber serialization is not
        supported; use :class:`AmberCV` instead."""
        if engine == "amber":
            raise NotImplementedError(
                "Direct Amber serialization of AngleCV is not supported. "
                "Wrap this CV with AmberCV to emit DISANG restraint files."
            )
        return {
            "cv_type": self.cv_type,
            "atom_1": self.atom_1_key,
            "atom_2": self.atom_2_key,
            "atom_3": self.atom_3_key,
            "window_target": window_target,
            "unit": self.unit,
        }

    def to_dict(self) -> Dict[str, Any]:
        return {
            "cv_type": self.cv_type,
            "name": self.name_,
            "atom_1_key": self.atom_1_key,
            "atom_2_key": self.atom_2_key,
            "atom_3_key": self.atom_3_key,
        }

    @classmethod
    def _from_dict_impl(cls, data: Dict[str, Any]) -> "AngleCV":
        return cls(
            atom_1_key=data["atom_1_key"],
            atom_2_key=data["atom_2_key"],
            atom_3_key=data["atom_3_key"],
            name=data.get("name", "angle_cv"),
        )


# ---------------------------------------------------------------------------
# DihedralCV
# ---------------------------------------------------------------------------

@_register_cv
class DihedralCV(CollectiveVariable):
    """Torsion (dihedral) angle defined by four atoms, measured in degrees.

    Attributes:
        name_:        Human-readable label.
        atom_1_key:   Key string for the first atom.
        atom_2_key:   Key string for the second atom.
        atom_3_key:   Key string for the third atom.
        atom_4_key:   Key string for the fourth atom.
    """

    cv_type: str = "dihedral"

    def __init__(
        self,
        atom_1_key: str,
        atom_2_key: str,
        atom_3_key: str,
        atom_4_key: str,
        name: str = "dihedral_cv",
    ) -> None:
        self.atom_1_key = atom_1_key
        self.atom_2_key = atom_2_key
        self.atom_3_key = atom_3_key
        self.atom_4_key = atom_4_key
        self.name_: str = name

    @property
    def name(self) -> str:
        return self.name_

    @property
    def unit(self) -> str:
        return "degree"

    def evaluate_on_structure(self, structure) -> float:
        """Compute the dihedral angle in degrees."""
        a1 = _resolve_atom(self.atom_1_key, structure)
        a2 = _resolve_atom(self.atom_2_key, structure)
        a3 = _resolve_atom(self.atom_3_key, structure)
        a4 = _resolve_atom(self.atom_4_key, structure)
        return float(a1.dihedral_with(a2, a3, a4))

    def serialize_for_engine(
        self,
        engine: str,
        work_dir: str,
        window_target: float,
    ) -> Dict[str, Any]:
        """Return metadata for *engine*.  Direct Amber serialization is not
        supported; use :class:`AmberCV` instead."""
        if engine == "amber":
            raise NotImplementedError(
                "Direct Amber serialization of DihedralCV is not supported. "
                "Wrap this CV with AmberCV to emit DISANG restraint files."
            )
        return {
            "cv_type": self.cv_type,
            "atom_1": self.atom_1_key,
            "atom_2": self.atom_2_key,
            "atom_3": self.atom_3_key,
            "atom_4": self.atom_4_key,
            "window_target": window_target,
            "unit": self.unit,
        }

    def to_dict(self) -> Dict[str, Any]:
        return {
            "cv_type": self.cv_type,
            "name": self.name_,
            "atom_1_key": self.atom_1_key,
            "atom_2_key": self.atom_2_key,
            "atom_3_key": self.atom_3_key,
            "atom_4_key": self.atom_4_key,
        }

    @classmethod
    def _from_dict_impl(cls, data: Dict[str, Any]) -> "DihedralCV":
        return cls(
            atom_1_key=data["atom_1_key"],
            atom_2_key=data["atom_2_key"],
            atom_3_key=data["atom_3_key"],
            atom_4_key=data["atom_4_key"],
            name=data.get("name", "dihedral_cv"),
        )


# ---------------------------------------------------------------------------
# AmberCV
# ---------------------------------------------------------------------------

def _eval_amber_expr(expr: Union[str, float], x: float) -> float:
    """Evaluate an Amber restraint expression that may contain the placeholder
    ``x`` (the window target value).

    Examples::

        _eval_amber_expr("x-1", 5.0)   → 4.0
        _eval_amber_expr("x+0.5", 5.0) → 5.5
        _eval_amber_expr(3.0, 5.0)     → 3.0  (plain float, unchanged)

    Args:
        expr: Either a plain numeric value or a string expression using ``x``.
        x:    The window target value substituted for the symbol ``x``.

    Returns:
        The evaluated float.
    """
    if isinstance(expr, (int, float)):
        return float(expr)
    if isinstance(expr, str) and "x" in expr:
        try:
            from sympy import sympify  # lazy import
            val = float(sympify(expr).subs("x", x))
            return round(val, 6)
        except Exception as exc:
            raise ValueError(
                f"Failed to evaluate Amber expression '{expr}' with x={x}: {exc}"
            ) from exc
    # Fallback: try plain float conversion
    return float(expr)


@_register_cv
class AmberCV(CollectiveVariable):
    """Collective Variable with first-class Amber DISANG/NMROPT restraint
    serialization.

    :class:`AmberCV` wraps an underlying :class:`CollectiveVariable` (typically
    :class:`DistanceCV`) and stores Amber-specific restraint parameters
    (``r1``–``r4``, ``rk2``, ``rk3``, ``ialtd``, …).  Calling
    :py:meth:`serialize_for_engine` with ``engine='amber'`` writes a DISANG
    ``.rs`` file into *work_dir* and returns the path.

    The ``r1``–``r4`` parameters accept either numeric values **or** string
    expressions containing ``x``, where ``x`` is replaced by *window_target*
    at serialization time.  For example ``r2='x'``, ``r1='x-1'`` produce a
    flat-bottom potential centred on the window target.

    Atom indices written to the DISANG file use :attr:`~enzy_htp.structure.Atom.idx`
    (the PDB serial number).  This is correct when the structure has been
    processed by tleap and the PDB serial numbers match the Amber topology
    atom ordering.

    Attributes:
        cv:            Wrapped :class:`CollectiveVariable` instance.
        name_:         Optional name override.  Falls back to ``cv.name``.
        amber_params_: Dict with Amber-specific restraint parameters.
        atom_keys_:    Atom key strings resolved at construction time.
                       For :class:`DistanceCV` this is ``[atom_1_key, atom_2_key]``.
    """

    cv_type: str = "amber"

    # Default Amber restraint parameters (flat-bottom style).
    DEFAULT_AMBER_PARAMS: Dict[str, Any] = {
        "r1": "x-1",
        "r2": "x",
        "r3": "x",
        "r4": "x+1",
        "rk2": 100.0,
        "rk3": 100.0,
        "ialtd": 0,
        "rs_filename": "cv.rs",
    }

    def __init__(
        self,
        cv: CollectiveVariable,
        amber_params: Optional[Dict[str, Any]] = None,
        name: Optional[str] = None,
    ) -> None:
        """Construct an :class:`AmberCV`.

        Args:
            cv:            The wrapped :class:`CollectiveVariable`.
            amber_params:  Dict overriding any :attr:`DEFAULT_AMBER_PARAMS` key.
                           Recognised keys:

                           - ``r1``, ``r2``, ``r3``, ``r4``: flat-bottom bounds
                             (numeric or ``'x±offset'`` expression).
                           - ``rk2``, ``rk3``: harmonic force constants
                             (kcal mol⁻¹ Å⁻²).
                           - ``ialtd``: Amber ``ialtd`` flag (0 or 1).
                           - ``rs_filename``: name of the DISANG file written
                             inside *work_dir*.
            name:          Optional human-readable label.  Defaults to
                           ``cv.name``.
        """
        self.cv = cv
        self.name_: str = name if name is not None else cv.name
        self.amber_params_: Dict[str, Any] = deepcopy(self.DEFAULT_AMBER_PARAMS)
        if amber_params is not None:
            self.amber_params_.update(amber_params)

    # --- CollectiveVariable interface -----------------------------------

    @property
    def name(self) -> str:
        return self.name_

    @property
    def unit(self) -> str:
        return self.cv.unit

    def evaluate_on_structure(self, structure) -> float:
        """Delegate to the wrapped CV."""
        return self.cv.evaluate_on_structure(structure)

    def evaluate_on_frame(self, frame) -> float:
        """Delegate to the wrapped CV."""
        return self.cv.evaluate_on_frame(frame)

    # --- Amber serialization --------------------------------------------

    def serialize_for_engine(
        self,
        engine: str,
        work_dir: str,
        window_target: float,
    ) -> Dict[str, Any]:
        """Write an Amber DISANG ``.rs`` file for the given *window_target*.

        The file is written to ``work_dir / rs_filename`` (where
        ``rs_filename`` comes from :attr:`amber_params_`).  For group-distance
        CVs (where the underlying ``cv`` does not have simple two-atom
        selections) you should subclass :class:`AmberCV` and override
        :py:meth:`_build_raw_rs_dict`.

        Args:
            engine:        Must be ``'amber'`` for DISANG serialization.
            work_dir:      Directory into which the ``.rs`` file is written
                           (created if it does not exist).
            window_target: The RC target value for this window (in ``cv.unit``).

        Returns:
            ``{'disang': <absolute-path-to-.rs-file>}``

        Raises:
            NotImplementedError: If *engine* is not ``'amber'``.
            RuntimeError:        If the underlying CV type is not supported.
        """
        if engine != "amber":
            raise NotImplementedError(
                f"AmberCV only supports engine='amber', got '{engine}'."
            )

        os.makedirs(work_dir, exist_ok=True)
        rs_filename: str = self.amber_params_.get("rs_filename", "cv.rs")
        rs_path = os.path.join(work_dir, rs_filename)

        raw_rs_dict = self._build_raw_rs_dict(window_target)
        self._write_disang_file([raw_rs_dict], rs_path)

        return {"disang": os.path.abspath(rs_path)}

    def _build_raw_rs_dict(self, window_target: float) -> Dict[str, Any]:
        """Build the raw restraint dict for a single ``&rst`` block.

        The method inspects the wrapped CV to obtain atom indices.  Currently
        :class:`DistanceCV`, :class:`AngleCV`, and :class:`DihedralCV` are
        supported.  For group-distance restraints (involving multiple atoms per
        group) call :meth:`_build_raw_rs_dict_from_atom_keys` directly.

        Subclasses may override this method for custom restraint layouts.

        Args:
            window_target: Window RC target; substituted for ``x`` in
                           ``r1``–``r4`` expressions.

        Returns:
            A dict suitable for :py:meth:`_write_disang_file`.
        """
        cv = self.cv
        # Collect ordered atom keys from the underlying CV
        if isinstance(cv, DistanceCV):
            atom_keys = [cv.atom_1_key, cv.atom_2_key]
        elif isinstance(cv, AngleCV):
            atom_keys = [cv.atom_1_key, cv.atom_2_key, cv.atom_3_key]
        elif isinstance(cv, DihedralCV):
            atom_keys = [cv.atom_1_key, cv.atom_2_key, cv.atom_3_key, cv.atom_4_key]
        else:
            raise RuntimeError(
                f"AmberCV._build_raw_rs_dict: unsupported wrapped CV type "
                f"'{type(cv).__name__}'.  Override _build_raw_rs_dict() in a "
                "custom AmberCV subclass for specialised restraints."
            )

        # Build the dict
        result: Dict[str, Any] = {}

        # iat: use pre-resolved integer indices if bind_topology() was called,
        # otherwise fall back to atom key strings (for documentation / debugging).
        resolved = self.amber_params_.get("_resolved_iat")
        result["iat"] = resolved if resolved is not None else atom_keys

        # r1–r4 and force constants
        for k in ("r1", "r2", "r3", "r4", "rk2", "rk3"):
            val = self.amber_params_.get(k)
            if val is not None:
                result[k] = _eval_amber_expr(val, window_target)

        # ialtd and other integer flags
        for k in ("ialtd", "ifvari", "ir6"):
            if k in self.amber_params_:
                result[k] = self.amber_params_[k]

        return result

    def bind_topology(self, structure) -> "AmberCV":
        """Resolve atom key strings to Amber (PDB serial) atom indices using a
        concrete :class:`~enzy_htp.structure.Structure`.

        Returns a new :class:`AmberCV` whose internal ``amber_params_`` stores
        resolved integer ``iat`` lists instead of key strings.  The original
        object is **not** modified.

        This method should be called before :py:meth:`serialize_for_engine`
        whenever the ``iat`` entries in the raw restraint dict should contain
        real integer indices rather than placeholder key strings.

        Args:
            structure: A :class:`~enzy_htp.structure.Structure` whose atoms
                       share the same serial numbering as the Amber prmtop
                       (i.e. the structure has been processed by tleap).

        Returns:
            A new :class:`AmberCV` with integer atom indices stored in
            ``amber_params_['_resolved_iat']``.
        """
        bound = deepcopy(self)
        raw = bound._build_raw_rs_dict(0.0)  # target doesn't matter here
        iat_keys = raw["iat"]  # list of key strings
        resolved = []
        for key in iat_keys:
            atom = _resolve_atom(key, structure)
            resolved.append(atom.idx)
        bound.amber_params_["_resolved_iat"] = resolved
        return bound

    # --- DISANG file writer (standalone, no interface dependency) -------

    @staticmethod
    def _write_disang_file(raw_rs_dict_list: List[Dict[str, Any]], rs_path: str) -> None:
        """Write a list of raw restraint dicts to a DISANG ``.rs`` file.

        The format follows section 27.1 of the Amber manual::

            &rst
             iat=1,100,
             r1=1.0, r2=2.5, r3=2.5, r4=4.0,
             rk2=100.0, rk3=100.0,
            &end

        Args:
            raw_rs_dict_list: List of restraint dicts (one per ``&rst`` block).
            rs_path:          Destination file path.
        """
        lines: List[str] = []
        for raw_dict in raw_rs_dict_list:
            lines.append("&rst")
            for k, v in raw_dict.items():
                if k.startswith("_"):
                    continue  # skip internal keys like _resolved_iat
                if k == "iat" or k.startswith("igr"):
                    # iat values: list of ints or strings
                    entries = ",".join(str(e) for e in v)
                    lines.append(f" {k}={entries},")
                elif isinstance(v, float):
                    lines.append(f" {k}={v:.4f},")
                else:
                    lines.append(f" {k}={v},")
            lines.append("&end")
        lines.append("")  # trailing newline

        with open(rs_path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines))

    # --- Structure integration ------------------------------------------

    def to_structure_constraint(
        self,
        structure,
        window_target: float,
        rs_filepath: str = "{mdstep_dir}/cv.rs",
    ):
        """Convert this :class:`AmberCV` to a :class:`StructureConstraint`
        compatible with the existing ``md_simulation`` infrastructure.

        The constraint's ``target_value`` is set to *window_target* so that
        the existing Amber mdin writer resolves ``r1='x-1'`` etc. correctly.

        Args:
            structure:      A :class:`~enzy_htp.structure.Structure` whose atoms
                            match the atom key strings in the wrapped CV.
            window_target:  RC target for this umbrella window (in ``cv.unit``).
            rs_filepath:    Path template for the DISANG file.  The placeholder
                            ``{mdstep_dir}`` is supported and resolved by the
                            existing Amber mdin writer.

        Returns:
            A :class:`~enzy_htp.structure.structure_constraint.StructureConstraint`
            instance with Amber parameters set appropriately.

        Raises:
            RuntimeError: If the wrapped CV type is not supported.
        """
        # Lazy import to avoid circular deps
        from enzy_htp.structure.structure_constraint import (
            DistanceConstraint,
            AngleConstraint,
            DihedralConstraint,
        )

        cv = self.cv
        # Build amber params with x-expressions (not yet evaluated)
        amber_params: Dict[str, Any] = {
            k: v
            for k, v in self.amber_params_.items()
            if k in ("r1", "r2", "r3", "r4", "rk2", "rk3", "ialtd", "ifvari", "ir6")
        }
        amber_params["rs_filepath"] = rs_filepath

        if isinstance(cv, DistanceCV):
            a1 = _resolve_atom(cv.atom_1_key, structure)
            a2 = _resolve_atom(cv.atom_2_key, structure)
            cons = DistanceConstraint(atoms=[a1, a2], target_value=window_target)
            cons.params["amber"].update(amber_params)
            return cons

        if isinstance(cv, AngleCV):
            a1 = _resolve_atom(cv.atom_1_key, structure)
            a2 = _resolve_atom(cv.atom_2_key, structure)
            a3 = _resolve_atom(cv.atom_3_key, structure)
            cons = AngleConstraint(atoms=[a1, a2, a3], target_value=window_target)
            cons.params["amber"].update(amber_params)
            return cons

        if isinstance(cv, DihedralCV):
            a1 = _resolve_atom(cv.atom_1_key, structure)
            a2 = _resolve_atom(cv.atom_2_key, structure)
            a3 = _resolve_atom(cv.atom_3_key, structure)
            a4 = _resolve_atom(cv.atom_4_key, structure)
            cons = DihedralConstraint(atoms=[a1, a2, a3, a4], target_value=window_target)
            cons.params["amber"].update(amber_params)
            return cons

        raise RuntimeError(
            f"AmberCV.to_structure_constraint: unsupported wrapped CV type "
            f"'{type(cv).__name__}'.  Supported: DistanceCV, AngleCV, DihedralCV."
        )

    # --- serialization --------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        return {
            "cv_type": self.cv_type,
            "name": self.name_,
            "wrapped_cv": self.cv.to_dict(),
            "amber_params": {
                k: v
                for k, v in self.amber_params_.items()
                if not k.startswith("_")
            },
        }

    @classmethod
    def _from_dict_impl(cls, data: Dict[str, Any]) -> "AmberCV":
        wrapped_cv = CollectiveVariable.from_dict(data["wrapped_cv"])
        return cls(
            cv=wrapped_cv,
            amber_params=data.get("amber_params"),
            name=data.get("name"),
        )


# ---------------------------------------------------------------------------
# PlumedCV
# ---------------------------------------------------------------------------

@_register_cv
class PlumedCV(CollectiveVariable):
    """Collective Variable declared via a raw PLUMED string.

    This class accepts either a complete PLUMED ``plumed.dat`` fragment or a
    single PLUMED action line.  It is intentionally minimal; more sophisticated
    parsing and conversion are left as future work.

    Attributes:
        plumed_string_: Raw PLUMED declaration(s).
        name_:          Human-readable label.
        unit_:          Unit string.
        label_:         PLUMED label used as the CV handle in the dat file.
    """

    cv_type: str = "plumed"

    def __init__(
        self,
        plumed_string: str,
        name: str = "plumed_cv",
        unit: str = "unknown",
        label: Optional[str] = None,
    ) -> None:
        """Construct a :class:`PlumedCV`.

        Args:
            plumed_string: Raw PLUMED action string or multi-line dat fragment.
                           Example: ``'d: DISTANCE ATOMS=1,10'``.
            name:          Human-readable label.
            unit:          Unit of the CV value (for documentation purposes).
            label:         PLUMED label used to reference this CV in
                           ``PRINT`` or ``METAD`` actions.  If *None* the
                           first word before ``:`` on the first line of
                           *plumed_string* is used (if present).
        """
        self.plumed_string_: str = plumed_string
        self.name_: str = name
        self.unit_: str = unit
        self.label_: Optional[str] = label or self._infer_label(plumed_string)

    @staticmethod
    def _infer_label(plumed_string: str) -> Optional[str]:
        """Try to extract the PLUMED label from the first line."""
        first_line = plumed_string.strip().splitlines()[0].strip()
        if ":" in first_line:
            return first_line.split(":")[0].strip()
        return None

    # --- CollectiveVariable interface -----------------------------------

    @property
    def name(self) -> str:
        return self.name_

    @property
    def unit(self) -> str:
        return self.unit_

    @property
    def plumed_string(self) -> str:
        """The raw PLUMED declaration string."""
        return self.plumed_string_

    @property
    def label(self) -> Optional[str]:
        """The PLUMED label for this CV."""
        return self.label_

    def evaluate_on_structure(self, structure) -> float:
        """Not implemented for PLUMED-only CVs.

        PLUMED CVs cannot be evaluated without running the PLUMED engine.
        Raises :py:exc:`NotImplementedError`.
        """
        raise NotImplementedError(
            "PlumedCV.evaluate_on_structure() is not implemented because "
            "PLUMED CVs require the PLUMED runtime to evaluate.  "
            "Compute the CV externally and pass the timeseries file to "
            "analysis functions."
        )

    def serialize_for_engine(
        self,
        engine: str,
        work_dir: str,
        window_target: float,
    ) -> Dict[str, Any]:
        """Write a ``plumed.dat`` fragment for *engine*.

        For ``engine='plumed'`` the PLUMED string is written to
        ``work_dir/plumed.dat`` and the path is returned.

        Args:
            engine:        Target engine; must be ``'plumed'``.
            work_dir:      Directory into which the dat fragment is written.
            window_target: Window target (appended as a comment for traceability).

        Returns:
            ``{'plumed_dat': <path>}``

        Raises:
            NotImplementedError: If *engine* is not ``'plumed'``.
        """
        if engine != "plumed":
            raise NotImplementedError(
                f"PlumedCV only supports engine='plumed', got '{engine}'."
            )

        os.makedirs(work_dir, exist_ok=True)
        dat_path = os.path.join(work_dir, "plumed.dat")

        content = (
            f"# window_target = {window_target}  unit = {self.unit}\n"
            + self.plumed_string_
            + "\n"
        )
        with open(dat_path, "w", encoding="utf-8") as fh:
            fh.write(content)

        return {"plumed_dat": os.path.abspath(dat_path)}

    # --- convenience constructors ---------------------------------------

    @classmethod
    def from_plumed(cls, plumed_string: str, **kwargs) -> "PlumedCV":
        """Construct a :class:`PlumedCV` from a raw PLUMED string.

        This is a thin alias for the constructor provided for API symmetry.

        Args:
            plumed_string: Raw PLUMED declaration(s).
            **kwargs:      Additional keyword arguments forwarded to
                           :py:meth:`__init__`.

        Returns:
            A new :class:`PlumedCV`.
        """
        return cls(plumed_string=plumed_string, **kwargs)

    # --- serialization --------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        return {
            "cv_type": self.cv_type,
            "name": self.name_,
            "unit": self.unit_,
            "label": self.label_,
            "plumed_string": self.plumed_string_,
        }

    @classmethod
    def _from_dict_impl(cls, data: Dict[str, Any]) -> "PlumedCV":
        return cls(
            plumed_string=data["plumed_string"],
            name=data.get("name", "plumed_cv"),
            unit=data.get("unit", "unknown"),
            label=data.get("label"),
        )
