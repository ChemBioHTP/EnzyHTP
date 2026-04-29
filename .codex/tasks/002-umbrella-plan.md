# Umbrella Sampling Implementation Plan for EnzyHTP

## Goal
Implement a modular umbrella sampling workflow that matches EnzyHTP architecture and coding conventions, with complete API coverage, external WHAM integration, robust testing, and machine-readable analysis outputs through dedicated typed output classes that can be converted to pandas DataFrames.

---

## 1) Data Model Layer: `CollectiveVariable`

### Objective
Introduce a reusable abstraction that captures the reaction coordinate definition, target value, and restraint parameters.

### Proposed location
- `enzy_htp/structure/collective_variable.py`
- Export in `enzy_htp/structure/__init__.py`

### Required capabilities
- Store:
  - `name: str`
  - `unit: str`
  - `target_value: float`
  - `force_constant: float`
  - optional engine-specific params
- Provide:
  - `build_constraint(topology: Structure) -> StructureConstraint`
  - `evaluate(structure_or_ensemble) -> float | np.ndarray`
- Preserve topology consistency checks aligned with existing constraint APIs.


- Idea: store code as strings?
- location of .so binary


- 1. run cv with md
- 2
- 

### Notes
- Keep this engine-agnostic.
- Reuse existing `StructureConstraint` creation patterns.

---

## 2) WHAM External Interface Layer

### Objective
Create a WHAM execution interface in `_interface`, consistent with EnzyHTP's interface design.

### Proposed files
- `enzy_htp/_interface/wham_interface.py` (new)
- Register in:
  - `enzy_htp/_interface/interface.py`
  - `enzy_htp/_interface/__init__.py`


### Required interface methods
- `set_parent(...)`
- executable/environment checks via config
- `write_metadata(...)`
- `run_wham(...)`
- `parse_pmf(...)`
- standardized error mapping/logging

**Important:** When generating WHAM metadata, ensure that the force constant value written is double the value used in AMBER (i.e., `force_constant_WHAM = 2 * force_constant_AMBER`). This is required because WHAM expects the force constant in a form that is double the AMBER restraint value. This nuance must be handled in the metadata writing logic and clearly documented in code and tests.

### Integration
- Use `eh_config["wham.*"]` defaults for executable, bins, tolerance, temperature, padding.
- Ensure environment checks include WHAM executable.

---

## 3) Umbrella Sampling Science API

### Objective
Provide a science API aligned with `equi_md_sampling` style while supporting multiple umbrella windows.

### Proposed location
- `enzy_htp/geometry/sampling.py` — add `umbrella_sampling()` alongside existing functions such as `equi_md_sampling()`


### Required improvements
- Prefer `CollectiveVariable` (or equivalent window-spec object) over raw callable-only input.
- Add explicit scheduling modes:
  - serial warm-start between windows
  - independent windows from same initial structure
- Return a structured result object (window metadata + ensembles + final structures + diagnostics), not a bare nested list/tuple.
- Keep tunable parameters explicit:
  - spring constant
  - window targets
  - production time
  - record period
  - temperature

**Note:** When passing force constants to the WHAM interface or generating metadata for WHAM, always apply the conversion: `force_constant_WHAM = 2 * force_constant_AMBER`. This ensures consistency with WHAM's expectations and avoids subtle analysis errors. This requirement should be reflected in both implementation and documentation.

### Output expectations
- Deterministic window directory layout.
- Machine-readable metadata for downstream analysis.

### Output class design: `UmbrellaSamplingResult`

The Science API should return a single `UmbrellaSamplingResult` object (not a bare list/tuple) that bundles all per-window data and global metadata. Its design should mirror the `MolDynResult` pattern: lazy file-backed access where trajectories are large, eager scalars elsewhere.

#### Proposed location
- `enzy_htp/geometry/umbrella_result.py` (new)
- Export via `enzy_htp/geometry/__init__.py`

#### Constituent class: `UmbrellaWindowResult`

Represents one biased window of the simulation.

| Field | Type | Description |
|-------|------|-------------|
| `window_id` | `int` | Zero-based index of this window |
| `collective_variable` | `CollectiveVariable` | CV definition (target value, force constant, unit, name) |
| `traj_file` | `str` | Path to trajectory file for this window |
| `traj_parser` | `Callable` | Lazy parser producing `List[Structure]` from `traj_file` |
| `log_file` | `str` | Path to MD log for this window |
| `log_parser` | `Callable` | Parser returning `Dict` of per-step metrics (energy, temperature, CV value) |
| `last_frame_file` | `str` | Path to final restart/frame (feed into next window or analysis) |
| `last_frame_parser` | `Callable` | Parser returning `Structure` |
| `work_dir` | `str` | Directory where window files reside |
| `metadata` | `Dict[str, Any]` | Catch-all for engine-specific extras (e.g., Amber restraint file path, `nstlim`, `dt`) |

Helper methods:
- `traj` → `List[Structure]` — calls `traj_parser(traj_file)` on demand
- `last_frame` → `Structure` — calls `last_frame_parser(last_frame_file)` on demand
- `log_data` → `Dict` — calls `log_parser(log_file)` on demand

#### Top-level class: `UmbrellaSamplingResult`

Aggregates all windows and global metadata.

| Field | Type | Description |
|-------|------|-------------|
| `windows` | `List[UmbrellaWindowResult]` | Ordered list of window results (index matches `window_id`) |
| `collective_variable_name` | `str` | Name of the shared CV axis (e.g., `"distance_Å"`) |
| `topology_file` | `str` | Shared topology/parameter file path |
| `topology_parser` | `Callable` | Parser returning `Structure` from topology |
| `scheduling_mode` | `str` | One of `"serial_warm_start"` or `"independent"` |
| `metadata` | `Dict[str, Any]` | Global metadata (engine, temperature, total windows, timestamp, etc.) |

Helper methods:
- `__len__()` — number of windows
- `__iter__()` — iterate over `UmbrellaWindowResult` objects
- `__getitem__(i)` — access window by index
- `window_targets` → `List[float]` — ordered target CV values
- `traj_files` → `List[str]` — ordered trajectory paths (convenience for WHAM metadata writing)
- `to_dataframe()` → `pd.DataFrame` — one row per window (columns: window_id, target_value, force_constant, unit, traj_file, last_frame_file, work_dir, plus flattened scalar metadata)
- `to_structure_ensemble(window_id: int)` → `StructureEnsemble` — wraps a single window's trajectory in a `StructureEnsemble` for geometry analysis

#### Design notes
- Immutable after construction; no setters. All mutable state lives in files on disk.
- `metadata` dicts use string keys only so they survive `json.dumps` trivially.
- `UmbrellaWindowResult` and `UmbrellaSamplingResult` should **not** import from `_interface` — they are pure data containers. Engine-specific result eggs (`AmberMDResultEgg`, etc.) are translated to `UmbrellaWindowResult` inside the Science API function.
- `to_dataframe()` column ordering must be deterministic (fixed list, not dict insertion order) so tests can compare without sorting.
- The `UmbrellaSamplingResult` is the natural input to the Section 4 analysis APIs (`wham_pmf`, `window_overlap_matrix`, etc.), which extract `traj_files` and `window_targets` from it.

#### Checklist additions
- [ ] `UmbrellaWindowResult` implemented with lazy accessors
- [ ] `UmbrellaSamplingResult` implemented with aggregation helpers and `to_dataframe()`
- [ ] Science API returns `UmbrellaSamplingResult` (not bare list)
- [ ] Unit tests: field access, lazy parsing, `to_structure_ensemble()`
- [ ] Analysis APIs accept `UmbrellaSamplingResult` directly (not raw file lists)

---

## 4) Umbrella Analysis API

### Objective
Complete analysis pipeline for density, PMF, and overlap diagnostics without plotting dependencies.

### Proposed location
- `enzy_htp/analysis/umbrella.py` (new)

### Scope boundary
- Do **not** add or require `matplotlib` in the core library implementation.
- Plot generation should remain external (e.g., user scripts/notebooks) and consume machine-readable outputs from APIs.
- Use a dedicated typed output class (or classes) for umbrella analysis results.
- Require explicit conversion helpers (e.g., `to_dataframe()`) so each output object can be converted to a pandas DataFrame for downstream tooling.

### Implementation requirement
- `wham_pmf(...)` must use `interface.wham` for all WHAM execution (no direct `subprocess` calls in analysis code).

### Required new diagnostics
- `window_overlap_matrix(...)`
- `detect_insufficient_overlap(...)`

### Required analysis outputs
- probability density machine-readable table/object
- PMF machine-readable table/object
- overlap diagnostic machine-readable report/object (warnings/threshold flags)

### Output representation conventions
- Public analysis APIs should return dedicated output classes with explicit fields/properties.
- Each output class must expose a deterministic `to_dataframe()` conversion path.
- Keep schema/field names stable and documented (e.g., coordinate/bin center, probability, free energy, count/weight, window id).
- Ensure deterministic ordering for reproducibility.
- Provide easy serialization paths via output-class helpers and/or DataFrame export (e.g., `result.to_dataframe().to_csv(...)`), while keeping plotting outside core APIs.

### Suggested output classes
- `UmbrellaDensityResult`
  - fields: `bin_centers`, `probabilities`, `counts`, `window_ids`, `metadata`
  - helper: `to_dataframe()`
- `UmbrellaPMFResult`
  - fields: `bin_centers`, `free_energies`, `uncertainties`, `counts`, `metadata`
  - helper: `to_dataframe()`
- `UmbrellaOverlapResult`
  - fields: `window_ids`, `overlap_matrix`, `threshold`, `flags`, `warnings`, `metadata`
  - helper: `to_dataframe()`

These classes should remain lightweight, typed containers with stable schemas that are easy to serialize and compare in tests.

---

## 5) Testing Strategy

### Unit tests
Add/expand tests for:
- `CollectiveVariable` construction, validation, constraint generation, and evaluation
- WHAM interface command assembly and parsing
- WHAM failure modes (missing executable, non-zero exit, malformed output)
- overlap diagnostics edge cases
- output-class schema/field validation, `to_dataframe()` conversion, and deterministic ordering

### Integration tests
Add/expand tests for:
- end-to-end umbrella workflow with analysis
- WHAM-enabled path (if executable available)
- fallback/skip behavior when WHAM unavailable
- reproduction of script behavior from:
  - `experiments/exp-umbrella/run18/template_md_only.py`
  - `experiments/exp-umbrella/run18/umbrella_figures.py`
- tests using the same data as `experiments/exp-umbrella/run18/template_md_only.py` which reproduce similar numerical results.

### Edge-case tests (required)
- insufficient phase-space overlap triggers diagnostic warning/error
- sparse bins and empty bins in WHAM setup
- non-monotonic or duplicated window targets
- constraint/topology mismatch

---

## 6) Completion Checklist (Measurable)

- [ ] `CollectiveVariable` implemented and exported
- [ ] WHAM interface implemented and registered in global interface
- [ ] umbrella science API updated to use CV-centered design
- [ ] density + PMF APIs finalized with dedicated output classes and `to_dataframe()` support
- [ ] overlap diagnostics implemented
- [ ] unit tests added for all new modules/features
- [ ] integration tests added for complete workflow
- [ ] regression test reproduces provided umbrella scripts
- [ ] all relevant tests pass in CI/local pytest

---

## 7) Suggested Rollout Order

1. Add `CollectiveVariable` type + unit tests
2. Add `WhamInterface` + registration + unit tests
3. Refactor analysis `wham_pmf` to use interface
4. Upgrade umbrella sampling API outputs and scheduling modes
5. Add overlap diagnostics reports
6. Add integration and regression tests
7. Run full targeted test suite and finalize docs/examples

---

## 8) Definition of Done

The feature is complete when:
1. Umbrella sampling can be launched from a reaction-coordinate abstraction with tunable restraint strength.
2. Probability density and PMF results are produced through public analysis APIs as dedicated output classes, each with DataFrame conversion support.
3. WHAM is integrated through `_interface` + `_config` patterns (no direct ad hoc execution in science/analysis APIs).
4. Overlap diagnostics identify insufficient sampling windows.
5. Unit + integration tests pass and cover all designed features and edge cases.
6. No `matplotlib` dependency is introduced by the umbrella sampling implementation.
7. Clean, from-scratch implementation follows EnzyHTP coding standards and architecture patterns.
