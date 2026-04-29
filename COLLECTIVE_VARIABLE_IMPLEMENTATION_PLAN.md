# Collective Variable (CV) Implementation Plan for EnzyHTP

## Goal
Provide a polymorphic, engine-agnostic Collective Variable (CV) abstraction and adapters so CVs can be declared programmatically or via PLUMED, serialized for MD engines (Amber, PLUMED-enabled engines), and consumed by analysis (probability densities, WHAM, plotting).

This design document mirrors the style of the repository's umbrella sampling plan and focuses specifically on the CV data model, engine adapters, serialization, analysis interoperability, and rollout strategy.

---

## 1) Data Model Layer: `CollectiveVariable` (CV)

### Objective
Introduce a reusable, typed abstraction that captures the reaction-coordinate definition, per-window target values, and restraint parameters while remaining engine-agnostic.

### Proposed location
- `enzy_htp/structure/collective_variable.py` (new)
- Export in `enzy_htp/structure/__init__.py`

### Core responsibilities
- Hold canonical CV metadata: `name`, `unit`, `selector`/`atom_indices` or expression, and engine-agnostic `parameters` (dict).
- Provide evaluation: compute a scalar RC for a `Structure` or a trajectory frame.
- Provide window generation: produce ordered `window_targets` from list or (start,end,step).
- Provide engine serialization hooks: produce engine-specific payloads (files or snippets) via `serialize_for_engine(engine, work_dir, window_target)`.
- Provide stable on-disk representation via `to_dict()` / `from_dict()`.

### Required API (sketch)
```py
class CVBase(ABC):
    name: str
    unit: str

    def generate_window_targets(self, start=None, end=None, step=None, list=None) -> List[float]: ...
    def serialize_for_engine(self, engine: str, work_dir: str, window_target: float) -> dict|str: ...
    def evaluate_on_structure(self, structure) -> float: ...
    def evaluate_on_frame(self, frame) -> float: ...
    def to_dict(self) -> dict: ...
    @classmethod
    def from_dict(cls, data: dict) -> 'CVBase': ...
```

### CV Subclasses (initial)
- `DistanceCV`: distance between two selections (residue- or atom-based). Units: Å.
- `AngleCV`: angle between three points. Units: degrees.
- `DihedralCV`: torsion angle between four points. Units: degrees.
- `AmberConstraintCollectiveVariable(CVBase)`: wraps Amber restraint generator (see section 2).
- `PlumedCollectiveVariable(CVBase)`: holds PLUMED declarations or converters.

### Acceptance criteria
- Instances produce deterministic numeric RC values for a given structure/frame.
- `to_dict()` must contain full information to reconstruct identical CV with `from_dict()`.

---

## 2) Amber Adapter

### Objective
Make Amber-friendly CVs by emitting the restraint (`.rst` / `DISANG` or inline md input) files and providing small md input snippets for job payload assembly.

### Proposed class
- `AmberConstraintCollectiveVariable(CVBase)` in `enzy_htp/structure/collective_variable.py`.

### Responsibilities
- Wrap existing `structure_constraint.create_group_distance_constraint`-style generators.
- Map engine-agnostic parameters to Amber restraint params (`r1..r4`, `rk2`, `rk3`, `ialtd`), exposing sensible defaults and unit checks.
- Implement `serialize_for_engine('amber', work_dir, window_target)` to write restraint files into `work_dir/{window_idx}/` and return a payload dict: `{'disang': path, 'md_snippet': '...'},` or path to restraint file.

### Window generation and parameterization
- `generate_window_targets` should return values in Å (canonical units). Conversion helper functions to engine units (if needed) are provided in the adapter.

### Integration hooks
- The umbrella MD launcher will look for `cv_payload = cv.serialize_for_engine('amber', mdstep_dir, target)` and include `cv_payload['disang']` in the MD input or copy into the run directory.

### Acceptance criteria
- The adapter writes valid restraint files matching patterns used elsewhere in the code (templating support for `{mdstep_dir}`), and the umbrella workflow can include produced files in MD jobs.

---

## 3) PLUMED Adapter & Declaration Support

### Objective
Support PLUMED-declared CVs as first-class objects: accept raw PLUMED strings or build PLUMED declarations from internal CV descriptions when possible.

### Proposed class
- `PlumedCollectiveVariable(CVBase)` in the same module or `enzy_htp/interface/plumed.py`.

### Responsibilities
- Accept either:
  - a PLUMED declaration string (e.g. `d: DISTANCE ATOMS=1,10`), or
  - an internal descriptor convertible to PLUMED (e.g., `DistanceCV` → `DISTANCE` directive).
- Implement `serialize_for_engine('plumed', work_dir, window_target)` to write `plumed.dat` fragments per window and return payload metadata.
- Provide `from_plumed(plumed_string)` best-effort parser for common CVs (document limitations).

### Edge cases
- Some PLUMED CVs have no simple mapping to internal descriptors (e.g., complex collective functions). For those, mark them as PLUMED-only and require PLUMED to be responsible for timeseries export; analysis functions then accept PLUMED timeseries files.

### Acceptance criteria
- PLUMED declarations can be embedded into MD jobs produced by `umbrella_sampling` and yield timeseries that analysis functions can consume.

---

## 4) Serialization & Analysis Interoperability

### Objective
Ensure analysis tools can compute probability densities and WHAM with respect to any `CVBase` instance.

### Requirements
- Stable, machine-readable CV serialization: JSON schema returned by `to_dict()` including engine-specific payloads and selectors.
- A helper `extract_reaction_coordinate` should accept either a `CVBase` or precomputed timeseries file; prefer `CVBase` and a trajectory parser to compute on-demand.
- Timeseries export function `export_timeseries_for_wham(timeseries, path)` to create WHAM-compatible files (frame index + RC value) and to attach necessary metadata (target distance, force constant).

### Implementation notes
- Keep unit canonicalization: densities and PMF are calculated in CV units (e.g., Å or degrees). Document conversions clearly.
- If CV is PLUMED-only and analysis requires RC values, require user to provide PLUMED timeseries file paths; implement `CVBase.from_plumed_timeseries(...)` minimal wrapper.

### Acceptance criteria
- `probability_density`, `wham_pmf`, and `plot_probability_density` accept `CVBase` or timeseries arrays and produce reproducible outputs.

---

## 5) Integration with Umbrella Sampling Workflow

### Objective
Wire `umbrella_sampling` (and `UmbrellaSamplingResult`) to accept any `CVBase` instance and use `serialize_for_engine` to prepare per-window files.

### Changes required
- `umbrella_sampling(..., cv: CVBase, ...)` — no API breaking change beyond accepting CVBase.
- For each window:
  1. call `payload = cv.serialize_for_engine(engine, mdstep_dir, target)`
  2. include `payload` files/snippets in the job directory
  3. ensure `UmbrellaWindowResult.metadata` includes `cv: cv.to_dict()` and `engine_payload: payload`

### Acceptance criteria
- Existing tests like `test/geometry/test_umbrella_integration.py` should run when using `AmberConstraintCollectiveVariable` with minimal changes.

---

## 6) Tests & Documentation

### Unit tests
- `test_structure/test_collective_variable.py`:
  - construction + `to_dict()`/`from_dict()` round-trip
  - `evaluate_on_structure` correctness for small canonical cases
  - `generate_window_targets` variants
  - `serialize_for_engine('amber')` output shape and basic content (file written)

### Integration tests
- Update `test/geometry/test_umbrella_integration.py` to exercise `AmberConstraintCollectiveVariable` and timeseries extraction.

### Documentation
- Add usage examples in `template/` and a short README section showing:
  - programmatic CV declaration (DistanceCV)
  - PLUMED declaration usage
  - exporting timeseries for WHAM

---

## 7) Minimum Implementation Sequence (first pass)
1. Implement `CVBase` and `DistanceCV` + `AmberConstraintCollectiveVariable`.
2. Add `generate_window_targets` and `serialize_for_engine('amber')` with file output.
3. Wire `umbrella_sampling` to accept `CVBase` and include payload in window metadata.
4. Add unit tests for CV behaviors and a small integration test using `test_umbrella_integration.py` snippet.
5. Document API and examples.

---

## 8) Completion Checklist
- [ ] `CVBase` implemented and exported
- [ ] `DistanceCV` + `AngleCV` + `DihedralCV` basic implementations
- [ ] `AmberConstraintCollectiveVariable` adapter implemented
- [ ] `PlumedCollectiveVariable` (skeleton) implemented
- [ ] `umbrella_sampling` accepts `CVBase` and stores payloads in `UmbrellaSamplingResult`
- [ ] Unit + integration tests added
- [ ] Documentation and usage examples added

---

## Definition of Done
The CV feature is complete when:
1. Users can declare CVs programmatically or via PLUMED and run umbrella sampling using the same high-level API.
2. Engine adapters produce the necessary files and metadata for MD job submission.
3. Analysis functions accept `CVBase` or precomputed timeseries and produce deterministic density/PMF outputs.
4. Tests and examples validate the main usage paths and the code follows EnzyHTP conventions.
