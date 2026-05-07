# Task: Segmented Local Velocity Reassignment and Scaling

Owner: `enzy_htp.geometry` + `enzy_htp._interface.amber_interface`  
Status: In progress. 

---

## Objective
Implement a segmented local-heating workflow in EnzyHTP that uses
**repeated local velocity reassignment / rescaling between short MD segments** to maintain
approximately sustained local heating in a selected region.

This task is explicitly **not** about true continuous local thermostatting during integration.
Amber should be treated as an external backend **used as-is**. Therefore the workflow must be
built entirely inside EnzyHTP by orchestrating short MD segments and modifying restart velocities
between segments.

The user-facing API should center on a function such as:

```python
md_energy_injection()
```

The implementation must fit the existing EnzyHTP framework instead of introducing a separate
workflow stack.

---

## Capability Boundary

This boundary must be explicit in code and documentation.

### What this task does support

- Standard Amber MD segments run through the existing EnzyHTP MD pipeline.
- Between segments, EnzyHTP reads restart coordinates/velocities, modifies the velocities of atoms
  in a selected driven region, writes a new restart, and launches the next segment.
- The supported perturbation modes should include:
  - **local heating by drawing new velocities from a Maxwell-Boltzmann distribution**
    for atoms in the driven region at a user-specified target temperature.
  - **velocity scaling** of the driven-region velocities to match a target kinetic temperature.

### What this task does not support

- True region-specific thermostatting during integration.
- Claims of exact local thermostat dynamics.
- Backend modifications to Amber.
- Silent presentation of this workflow as physically equivalent to continuous local thermostatting.

The workflow should be documented as a **segmented local-heating protocol** that may be
useful for pathway exploration, qualitative transport studies, and controlled local perturbation.

---

## Scientific Intent

The intended physical picture is:

- equilibrate the system conventionally,
- run short production MD segments,
- after each segment, reheat a chosen driven region by reassigning or rescaling velocities,
- continue this cycle over many segments.

This yields a driven trajectory in a **heuristic** sense over many small
individual production MD segments.

This workflow is appropriate for:

- exploratory vibrational energy transport studies,
- qualitative pathway probing,
- testing response of observables under repeated local heating,
- generating driven trajectories for downstream analysis.

This workflow is not appropriate when the project requires:

- exact local thermostat dynamics,
- formal claims about the stationary ensemble of a continuously driven stochastic process,
- strict equivalence to a patched MD integrator.

---

## Existing EnzyHTP Capabilities To Build On

The design should reuse what already exists in `enzy_htp/`.

### Existing MD workflow contracts

- `enzy_htp.geometry.sampling.md_simulation()` is the engine-agnostic MD orchestration API.
- `MolDynParameterizer`, `MolDynParameter`, `MolDynStep`, and `MolDynResult` already define the
  existing MD execution contracts.
- `equi_md_sampling()` already embodies the conventional equilibration path.
- `deployable_md_simulation()` and `deployable_equi_md_sampling()` already establish the pattern
  for submission-ready workflow assembly.

### Existing Amber extension seam

- `AmberInterface.build_md_parameterizer()`
- `AmberInterface.build_md_step()`
- restart / trajectory parsers already used by Amber MD result handling
- Amber mask generation via `AmberInterface.get_amber_mask()`

### Existing structure and selection tools

- `Structure`, `Atom`, `Residue`, `StruSelection`
- `select_stru()` for pymol-like user-facing selection strings

The new workflow should normalize all region-like inputs to `StruSelection` before any perturbation
logic is applied.

### External reference and citation requirements

Any new numerical constant, unit conversion, or Amber-specific restart/MD control behavior added for
this task must be documented with at least one source.

Required standard:

- general physical constants and conversions should cite a reputable physical reference,
- atomic-mass usage should cite a clear source for the mass convention being used,
- Amber-specific restart and mdin behavior should cite the Amber manual and, when possible, the
  specific relevant section.

These citations may live in code comments, parser docstrings, or adjacent developer
documentation, but they must be easy for a reviewer to trace.

---

## Design Requirements

### 1. Extend the current MD stack, do not fork it

The new workflow should be built on:

- `geometry.md_simulation()`
- standard Amber MD steps
- ordinary Amber restart files

The new code should assemble a sequence of MD segments and inter-segment perturbation steps, not a
parallel execution framework.

`md_simulation()` should be treated as the architectural reference and, where practical, the
general template for submission handling, run orchestration, and result packaging. The
energy-injection implementation should not introduce a second parallel orchestration stack with
separate child/shrapnel-style job logic when the same responsibilities can be expressed by
reusing or extending the current MD workflow pattern.

### 1a. Submit one segmented MD worker job per replica, not per segment

This workflow must not submit each production segment as a separate cluster job.

For realistic segment counts, repeated queueing would dominate walltime and make the workflow
operationally impractical. The EnzyHTP design should therefore keep the segment loop inside one
submitted worker job for each replica:

- the manager process or manager Slurm job prepares the workflow and submits MD worker jobs,
- the manager itself should not run GPU MD segments locally,
- each submitted worker job requests the GPU resources it needs,
- inside that worker job, all production segments for that replica run sequentially,
- restart velocities are mutated between segments inside the same running worker job.

The unit of GPU cluster submission should therefore be one whole segmented-MD replica workflow, not
an individual production segment. A higher-level manager job may still be used, but only to submit
and monitor these replica worker jobs.

### 2. Keep user semantics backend-independent where possible

Expose concepts such as:

- driven region
- target driving temperature
- perturbation interval
- per-segment MD length
- number of production cycles

Do not expose Amber implementation details in the public API unless required.

In particular:

- `cluster_job_config` should be handled as a config object or dictionary, not as a user-visible
  `"default"` sentinel within this workflow contract,
- scheduler-specific fields such as `res_keywords` must not become part of the science API surface,
- resource settings must remain user-configurable through the normal configuration objects rather
  than being buried inside internal helper functions.

### 3. Use current EnzyHTP region types

Accepted region inputs should include:

- `StruSelection`
- selection strings
- residue object(s)
- atom object(s)

Normalize these to `StruSelection` first.

### 4. Preserve current packaging style

Preferred placement:

```text
enzy_htp/
├── geometry/
│   └── sampling.py
└── _interface/
    └── amber_interface.py
```

Do not create a large new subpackage unless implementation size justifies it.

---

## Recommended Public API

```python
def md_energy_injection(
    stru,
    param_method=None,
    engine="amber",
    work_dir="./MD_SEGMENTED_INJECTION",
    prod_time=10.0,
    segment_time=0.05,
    temperature=300.0,
    driven_region=None,
    injection_mode="maxwell_reassign",
    driving_temperature=None,
    include_hydrogens=False,
    parallel_runs=1,
    parallel_method=None,
    cluster_job_config=None,
    record_period=None,
    return_observables=False,
    **kwargs,
):
    ...
```

```python
def deployable_md_energy_injection(
    stru,
    param_method=None,
    engine="amber",
    work_dir="./MD_SEGMENTED_INJECTION",
    cluster_job_config=None,
    **kwargs,
):
    ...
```

### API notes

- `engine` should support only `"amber"` in v1.
- `driven_region` is required.
- `driving_temperature` is required.
- `segment_time` is the length of each production segment between perturbations.
- `prod_time` is the total production time, which determines the number of segments.
- `injection_mode` should support both:
  - `"maxwell_reassign"`
  - `"velocity_scale"`
- `injection_mode="maxwell_reassign"` should be the default.
- `md_energy_injection(..., parallel_method=None)` should perform immediate execution in the
  current process.
- deployable cluster execution should be exposed through `deployable_md_energy_injection()`.
- if a convenience submission wrapper is later added, it should still submit only one cluster job
  per replica workflow.
- `cluster_job_config` should accept the normal EnzyHTP config forms (config object or dictionary)
  and should not rely on a string sentinel such as `"default"` in this workflow API.
- user-facing API signatures should not expose scheduler-detail knobs such as `res_keywords`
  directly; those belong inside configuration objects.

### Naming guidance

Prefer names like:

- `md_energy_injection`

Avoid names that imply exact local thermostatting.

---

## Recommended Internal Flow

The workflow should look like:

1. validate engine,
2. build or accept `param_method`,
3. normalize the driven region to `StruSelection`,
4. run standard equilibration,
5. split production into many short MD segments,
6. after each segment:
   - read the segment restart,
   - modify driven-region velocities,
   - zero net momentum if appropriate,
   - write a new restart,
7. launch the next MD segment from that modified restart,
8. collect segment outputs into a coherent result object or result list.

The implementation should reuse ordinary Amber MD steps for the production segments rather than
inventing a nonstandard execution model.

### Recommended cluster execution model

Cluster execution should still submit one whole replica workflow per allocation, but the submission
and orchestration pattern should be aligned with `md_simulation()` rather than built as a separate
child/shrapnel framework.

Required behavior:

1. `deployable_md_energy_injection()` must prepare a submission-ready workflow using the same
   overall handling style as the standard MD stack.
2. The segmented production loop must run inside one allocated worker execution context per replica.
3. No dedicated shrapnel-style orchestration layer, redundant child-main stack, or per-segment
   scheduler submission chain should be introduced for this task.
4. Any workflow-specific wrapper that remains should be minimal, easy to read, and clearly justified
   by the need to mutate restart velocities between completed MD segments.

The main design goal is not merely to avoid per-segment queueing, but also to avoid duplicating the
submission/orchestration logic that already exists elsewhere in `geometry.sampling`.

### Why `md_simulation(..., parallel_method="cluster_job")` is not enough by itself

The existing `cluster_job` path is appropriate when all steps can be declared as cluster jobs up
front. That does not match segmented energy injection, because:

- each next production segment depends on Python-side restart editing after the previous segment,
- the next runnable command cannot be finalized until the previous segment has completed,
- queueing each segment independently would create unacceptable scheduler overhead.

The segmented workflow should therefore not be implemented as a long sequence of ordinary EnzyHTP
MD job submissions. However, that limitation is not a justification for creating a fully separate
submission framework. The implementation should instead extend or adapt the existing
`md_simulation()` orchestration pattern as narrowly as possible.

---

## Velocity Reassignment Strategy

### Supported mode 1: Maxwell-Boltzmann redraw

For each perturbation cycle:

- identify atoms in the driven region,
- sample new Cartesian velocities from a Maxwell-Boltzmann distribution consistent with the target
  driving temperature,
- assign those velocities to the driven atoms,
- optionally remove driven-region center-of-mass drift or whole-system net momentum afterward.

This should be the default mode because it is:

- simple,
- explicit,
- easy to document honestly,
- less ad hoc than deterministic one-shot scaling.

### Supported mode 2: deterministic velocity scaling

For each perturbation cycle:

- compute current kinetic temperature of the driven region,
- rescale velocities to match the target temperature.

This mode should be supported because it is:

- simple,
- inexpensive,
- deterministic,
- useful for controlled perturbation studies.

It should not be the default if Maxwell redraw is available.

### Mode semantics

The public API should make the distinction explicit:

- `injection_mode="maxwell_reassign"`:
  redraw driven-region velocities from a Maxwell-Boltzmann distribution
- `injection_mode="velocity_scale"`:
  preserve current velocity directions but rescale magnitudes to match the target kinetic temperature

Do not silently switch between these modes.

### Recommended MD workflow for `velocity_scale`

The `velocity_scale` workflow should resemble `equi_md_sampling()` in overall shape, but the
production stage must be decomposed into an ordered sequence of short restart-based MD runs.

Recommended sequence:

1. **Minimization**
   - one minimization step, matching the current EnzyHTP convention of relaxing the prepared
     structure before dynamics.
2. **Heating equilibration**
   - one short NVT heating step that brings the system from low temperature to the target bulk
     temperature.
3. **Density / pressure equilibration**
   - one or two short NPT equilibration steps, following the current `equi_md_sampling()` pattern.
4. **Segmented production**
   - replace the single long production run with `n_segments = ceil(prod_time / segment_time)`
     consecutive production steps.
   - segment 0 starts from the final equilibration restart.
   - after segment `i` finishes, EnzyHTP must read the produced restart, extract coordinates and
     velocities, scale the driven-region velocities to the requested `driving_temperature`, write a
     new restart, and use that restart as the input for segment `i + 1`.
   - every production segment after the first must therefore run with restart semantics enabled and
     must consume velocities from the previous segment's modified restart.
   - when running on a cluster, all of these segment runs should execute inside the same submitted
     GPU job.

This is intentionally similar to `equi_md_sampling()` in that it still uses:

- one parameterization,
- one minimization stage,
- conventional equilibration before production,
- ordinary Amber MD steps produced by `build_md_step()`.

It differs from `equi_md_sampling()` in two required ways:

- production is a chain of short consecutive MD segments rather than one single long production
  step,
- production restarts must preserve, output, and re-read velocities at each segment boundary.

For Amber-facing step construction, the production-step requirements should be explicit:

- segment steps must write restart files that contain both coordinates and velocities,
- segment 0 may start from the equilibration restart, but segment `i > 0` must start from the
  modified restart generated after segment `i - 1`,
- each segment should still emit its normal `.out`, `.nc`, and restart artifacts so downstream
  EnzyHTP parsing remains unchanged.
- cluster deployment must keep these segment executions local to the workflow driver rather than
  resubmitting them through the scheduler.

---

## Restart Handling Requirements

This task depends on manipulating restart velocity data.

The implementation therefore needs:

- a reliable Amber restart parser that can read coordinates and velocities,
- a reliable Amber restart writer built around the same parsed representation,
- consistency with existing EnzyHTP Amber restart parsing utilities.

The preferred design is to introduce an explicit Amber restart parser abstraction rather than
embedding ad hoc `read_from_rst()` / `write_to_rst()` logic directly on the interface.

Implications:

- an `AmberRestart` data container should not remain as a separately maintained special-case API if
  the parser object can own the parsed representation cleanly,
- restart parsing and writing should become parser responsibilities,
- restart mutation logic for energy injection should operate on the parser-level representation.

The perturbation logic should operate on restart data, not on trajectory-only snapshots.

---

## Checkpointing and Resubmission Requirements

Because the workflow may run for a long time inside one GPU allocation, EnzyHTP should make it
restartable.

The workflow driver should checkpoint at least:

- completed equilibration status,
- completed production segment index,
- current restart file path,
- perturbation metadata already applied,
- any per-segment result metadata needed for final packaging.

If the allocated job ends before the workflow completes, a rerun of the driver should resume from
the last completed segment rather than restart the entire calculation from scratch.

---

## Numerical Safeguards

The implementation should include:

- exclusion of hydrogens by default,
- optional center-of-mass momentum removal after perturbation,
- checks that driven region is non-empty,
- validation that `segment_time > 0`,
- validation that `prod_time >= segment_time`,
- validation that the number of segments is at least 1.

The workflow should also warn if:

- segment lengths are extremely short,
- driving temperature is extremely high,
- the number of segments is extremely large.

The corresponding constants and conversion factors used in these safeguards and in the restart
velocity perturbation code must be documented with citations as noted above.

---

## Result Packaging

Prefer reusing the current `md_simulation()` return contract as much as possible.

If a wrapper is needed, it should add only:

- perturbation protocol metadata,
- per-segment result grouping,
- optional derived observables.

Do not duplicate trajectory/restart paths that already exist in `MolDynResult` unless there is a
clear benefit.

The implementation should also prioritize readability:

- keep the geometry-level orchestration legible to a human reviewer,
- keep workflow state objects simple and well named,
- avoid deep, redundant layers of submission wrappers,
- prefer small helpers with explicit responsibilities over large mixed-purpose functions.

---

## Analysis Guidance

The first implementation should stay modest.

Good initial follow-up quantities:

- segment index vs perturbation metadata,
- driven-region RMSD / structural response,
- trajectory-based downstream observables already supported by EnzyHTP,
- optional energy-injection bookkeeping estimated from velocity changes.

Avoid overpromising:

- exact local temperature measurements,
- formal thermodynamic interpretation beyond the documented heuristic protocol,
- mdout-based thermodynamic parsing that EnzyHTP does not yet support.

---

## Testing Requirements

Tests should be split across:

- `test/geometry/`
- `test/_interface/`

The definition of done is not just object construction. The implementation is only complete when
the workflow is shown to run and produce correct artifacts.

### 1. First-priority API tests: `md_energy_injection()` behavior must be exercised directly

The first and most important tests for this task are high-level pytest tests of
`md_energy_injection()` itself and its user-facing feature combinations.

These tests should be modeled after the manual validation that was done in the scripts under:

- `/home/larkinsd/workshop_new/Energy_injection_tests/energy_injection_001/`
- `/home/larkinsd/workshop_new/Energy_injection_tests/energy_injection_002/`
- `/home/larkinsd/workshop_new/Energy_injection_tests/energy_injection_003/`
- `/home/larkinsd/workshop_new/Energy_injection_tests/energy_injection_004/`

but expressed in a proper pytest-enabled way.

Required API-level feature coverage should include combinations of:

- `injection_mode="maxwell_reassign"`
- `injection_mode="velocity_scale"`
- multiple `driven_region` input forms and selections
- different `driving_temperature` values
- hydrogen-included vs default-hydrogen-excluded behavior when relevant

The main assertion strategy should be behavioral rather than purely structural:

- for every injection step, the kinetic energy of the selected driven region after perturbation
  should be significantly larger than before perturbation.

It is recommended to define one reusable helper for these tests that:

- identifies the driven-region atoms,
- reads the pre-injection and post-injection restart velocities,
- computes the selected-region kinetic energy,
- asserts that the post-injection value exceeds the pre-injection value by a meaningful margin.

This helper should then be reused across the main `md_energy_injection()` test matrix.

### 2. Regression test: standard MD must still work

The first high-level test gate is regression safety for the existing EnzyHTP framework.

At least one existing ordinary Amber MD workflow must still complete successfully through the
normal EnzyHTP path:

- `MolDynParameterizer -> MolDynStep -> geometry.md_simulation()`

This regression coverage should verify:

- parameterization still succeeds,
- standard MD input files are generated correctly,
- the MD job still runs successfully,
- expected output artifacts are still produced.

This is the first required done criterion. If normal MD is broken by the energy-injection changes,
the task is not complete.

### 2a. Regression test: deployable single-job execution must work

The energy-injection workflow must also be validated in its intended cluster execution mode.

At least one deployable workflow assembly test should verify that:

- one submission-ready cluster job is produced for one workflow replica,
- the generated job runs a Python workflow driver rather than one scheduler submission per segment,
- segment execution is represented inside the driver logic, not as a long scheduler-side job list.

### 3. End-to-end energy-injection MD must run

The second high-level test gate is an end-to-end segmented energy-injection workflow run, not only
mocked assembly.

At least one `md_energy_injection()` workflow must complete successfully for:

- `injection_mode="maxwell_reassign"`

and, if both modes are claimed in the implementation, also for:

- `injection_mode="velocity_scale"`

This coverage should verify:

- equilibration runs before the segmented driven production,
- production segments are executed sequentially,
- restart velocities are modified between segments,
- the segmented production is executed under one workflow-level cluster allocation in deployable
  mode,
- the workflow completes and produces intended output files.

### 4. Input and output artifact validation

Tests must validate the actual files produced by both standard MD and energy-injection MD.

Required checks should include:

- MD input files are written and contain the expected settings for the relevant step,
- output files such as `.out`, `.nc`, and restart files are produced where expected,
- these files are non-empty and large enough to indicate a real calculation rather than an empty or
  stub output,
- output files contain the minimum information needed for downstream EnzyHTP usage,
- trajectories and restart files remain parseable by the existing EnzyHTP machinery.

Do not treat file existence alone as sufficient evidence that the workflow works.

### 5. Region selection and perturbation-behavior tests

Tests must focus on whether the workflow behaves as intended, not only whether values are stored in
objects.

Required coverage should include:

- user-facing structure selection resolves correctly to `StruSelection`,
- hydrogen filtering defaults behave as intended,
- empty selections fail clearly,
- driven-region-only setups are tested,
- velocity modification is applied only to the intended selected region,
- unselected atoms are not modified except for any documented whole-system momentum correction.

When both injection modes are implemented, these checks should be performed for:

- `maxwell_reassign`
- `velocity_scale`

### 6. Restart mutation integrity tests

Because the workflow depends on modifying restart velocities, tests must verify parser-based restart
integrity.

Required checks should include:

- restart files can be parsed before perturbation,
- parsed restart data can be written back out after perturbation,
- modified restart files can be parsed again successfully,
- coordinates are preserved when only velocities are intended to change,
- velocity arrays change only in the selected region, subject to any documented momentum-removal
  step.

### 7. Interface-level and geometry-level split

Focus interface-level tests on `test/_interface/`.

Add tests for:

- region normalization,
- velocity perturbation helper logic for both supported modes,
- restart parser read/modify/write integrity,
- hydrogen filtering default behavior.

Focus geometry-level tests on `test/geometry/`.

Add tests for:

- `md_energy_injection()` argument validation,
- correct assembly of equilibration + segmented production steps,
- correct invocation of perturbation between segments,
- correct submission/orchestration behavior without a redundant child/shrapnel layer,
- correct single-job deployable assembly,
- end-to-end production of the expected workflow artifacts.

### 8. Amber mdin safety tests

Because `ntxo=1` is dangerous as a global change, tests must explicitly verify that restart-output
format changes are gated behind an energy-injection-specific handle.

Required checks should include:

- `AmberInterface.build_md_step()` / `AmberMDStep.md_config_dict` can represent whether a step is
  part of the energy-injection workflow,
- when `ascii_rst=True`, the generated mdin includes `ntxo=1`,
- when `ascii_rst=False` or unspecified, `ntxo=1` is absent,
- ordinary non-energy-injection MD workflows retain their previous restart-output behavior.

This is a required regression guard, not an optional convenience test.

---

## Recommended Implementation Order

1. revise the task architecture so the segmented workflow follows the existing
   `md_simulation()` handling pattern instead of introducing redundant child/shrapnel logic,
2. add an Amber restart parser/writer abstraction for restart read/modify/write operations,
3. add or clarify citations for physical constants, masses, and Amber-specific behavior,
4. implement region normalization and validation helpers if missing,
5. implement a velocity-perturbation helper for `maxwell_reassign`,
6. implement a velocity-perturbation helper for `velocity_scale`,
7. add an explicit md-step configuration handle such as `ascii_rst` so `ntxo=1` is opt-in for
   energy-injection-compatible steps only,
8. implement the geometry-level orchestration function for segmented production with restart-based
   velocity carryover,
9. add parser-level restart mutation tests and mdin safety tests,
10. add high-level pytest API tests modeled after the existing manual energy-injection examples,
11. add end-to-end workflow tests on a small system.

---

## Definition of Done

This task is done when:

1. EnzyHTP can still perform at least one full ordinary Amber MD workflow without regression.
2. `md_energy_injection()` can perform equilibration followed by segmented driven production.
3. Driven-region velocities are modified between segments using the documented supported injection
   mode(s).
4. The primary pytest coverage directly exercises `md_energy_injection()` across both perturbation
   modes, multiple region specifications, and multiple driving temperatures, with a reusable helper
   that verifies selected-region kinetic energy increases after each injection step.
5. Restart parsing and writing for this workflow are handled through a clear Amber restart parser
   abstraction rather than ad hoc interface-side read/write helpers.
6. `ntxo=1` is enabled only for energy-injection-compatible restart-writing steps through an
   explicit step/config handle such as `ascii_rst`, and ordinary MD workflows are protected by
   regression tests.
7. `cluster_job_config` handling remains type-consistent as a config object or dictionary, not a
   user-visible `"default"` sentinel in this workflow contract.
8. User-facing resource configuration remains exposed through configuration objects rather than
   leaking `res_keywords`-style scheduler details into the science API surface.
9. Produced input and output artifacts (`.out`, `.nc`, restart files, and relevant input files) are
   written as intended, are non-empty, and remain parseable by EnzyHTP.
10. Structure selection works correctly, including hydrogen filtering and empty-selection rejection.
11. The perturbation acts selectively on the intended driven region.
12. The workflow is clearly documented as segmented heuristic energy injection rather than
    continuous local thermostatting.
13. Cluster deployment submits the workflow as one GPU job per replica and executes segment chaining
    inside that allocation rather than queueing one job per segment.
14. The geometry-level and interface-level implementation are readable enough that a reviewer can
    trace the workflow without stepping through a redundant orchestration stack.
15. Focused interface-level and geometry-level tests pass.
