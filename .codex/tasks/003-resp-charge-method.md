# 003 - RESP Charge Method Support
# IMPORTANT: remind me to make sure of "definition of done" for this milestone.

## Scope

This task tracks implementation of `charge_method="RESP"` in the Amber `MolDynParameterizer` path. It is separate from MCPB.py metal-center automation, but the MCPB completion test depends on it because MCPB step `3b` uses RESP charges for renamed metal-center residues.

The target is not to redesign every charge workflow. The target is to make RESP a real, explicit charge method instead of a partially supported option or an accidental fallback to AM1-BCC behavior.

## Required Behavior

When a user builds an Amber parameterizer with:

```python
param_worker = interface.amber.build_md_parameterizer(
    ...,
    charge_method="RESP",
)
params = param_worker.run(stru)
```

the ligand, modified-residue, and MCPB-relevant charge workflow must use RESP-derived charges. It must not silently use AM1-BCC or another charge method.

## Design Constraints

- `MolDynParameterizer` remains structure-independent until `run(stru)`.
- Structure-specific charge/spin data must be read from the `Structure` passed to `run(stru)`, including values assigned through `assign_ncaa_chargespin()`.
- RESP artifacts should follow existing parameterizer temporary-directory and NCAA-library conventions where possible.
- Missing external QM/RESP outputs should fail early with an actionable error unless the implementation is explicitly configured to generate them.

## Development Tasks

### 1. Audit Existing Charge Method Routing

Identify where `charge_method` is interpreted for ligand and modified-residue parameterization.

Done when:

- `RESP` and `AM1BCC` take distinct code paths.
- Unsupported charge methods raise clear errors before external tools are invoked.
- No MCPB workflow can silently fall back from `RESP` to `AM1BCC`.

### 2. Implement RESP Artifact Handling

Teach the Amber parameterization flow how to provide or locate RESP charge artifacts for noncanonical residues and MCPB models.

Responsibilities:

- Define the expected RESP input/output files and where they live.
- Connect RESP charges to generated mol2 files for ligands and modified residues.
- Preserve charge/spin assignments from the input `Structure`.
- Report missing or inconsistent RESP artifacts with file paths and residue/model identifiers.

Done when:

- `charge_method="RESP"` produces parameter files with RESP charges for supported fixtures.
- Missing RESP artifacts fail with actionable diagnostics.

### 3. Support MCPB Step `3b` Needs

The MCPB task uses `mcpb_step3_method="3b"`, which fits RESP charges for the metal-center model.

Responsibilities:

- Provide MCPB step 3 with the log/charge artifacts it expects.
- Derive small/large model charge and spin from structure-level charge/spin data where possible.
- Keep MCPB-specific structure derivation inside the MCPB task; this task should only supply the RESP charge-method behavior it needs.

Done when:

- The MCPB 3PZW completion test can use `charge_method="RESP"` without charge-method fallback.
- RESP-related failures happen before tleap and point to the missing MCPB/RESP artifact.

## Completion Criteria

This task is complete when all of the following are true:

1. `charge_method="RESP"` is an implemented, explicit Amber parameterizer path.
2. `charge_method="RESP"` never silently falls back to AM1-BCC.
3. Required RESP artifacts are either generated or located through documented parameterizer/library conventions.
4. Missing RESP artifacts produce clear EnzyHTP errors.
5. The MCPB 3PZW completion test can depend on RESP behavior without owning RESP implementation details.

## Explicit Non-Goals

- Full automation of Gaussian/GAMESS job submission on HPC.
- Redesigning the entire NCAA parameterization API.
- Implementing MCPB.py metal-center discovery or tleap integration; that remains in `002-metal-support.md`.
