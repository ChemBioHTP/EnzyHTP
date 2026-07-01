# 002 - Metal Center Support Through MCPB.py

## Scope

This task focuses on full Amber bonded metal-center parameterization through `MCPB.py`.

The current milestone is not about general PDB parsing, chain handling, preparation, mutation, or the already-supported manual/nonbonded metal paths. Those layers may receive small integration hooks only when MCPB.py output needs to flow through the Amber parameterizer.

Current status in the codebase:

- `MetalUnit` and `Structure.metalcenters` already exist.
- `AmberParameterizer` already supports user-supplied `additional_tleap_lines`.
- `AmberParameterizer._parameterize_metalcenter()` is currently a TODO.
- `metal_involved_main.py` demonstrates a manually assembled bonded model by loading ion parameters and manually adding `bond` commands.

Target status:

- A user can provide a `Structure` containing a metal center and enough MCPB configuration, then `AmberParameterizer.run()` automatically drives MCPB.py, loads the generated MCPB mol2/frcmod/PDB/tleap inputs, and produces Amber `prmtop`/`inpcrd` equivalent to a manually run MCPB.py workflow.

## MCPB.py Workflow To Support

The local supported environment is:

```bash
source ~/bin/miniconda3/bin/activate new_enzy_htp
```

The available MCPB.py is AmberTools 23.6, `MCPB.py` version 7.0.

MCPB.py takes an input file and a step number:

```bash
MCPB.py -i input_file -s step_number [--logf gaussian_or_gamess_log] [--fchk gaussian_fchk]
```

Required input variables:

- `ion_ids`: atom IDs of metal ions in the original PDB.
- `ion_mol2files`: mol2 files for the metal ions.
- `original_pdb`: original PDB used by MCPB.py.

Important optional variables:

- `ion_info`: required for step `4n2`; residue name, atom name, element, charge.
- `add_bonded_pairs`: explicit bonded atom pairs, useful when auto-detection is insufficient.
- `additional_resids`: extra residue IDs to include in the metal model.
- `cut_off`: coordination cutoff, default 2.8 A in MCPB.py.
- `force_field`: Amber force field, default `ff19SB` in MCPB.py.
- `frcmod_files`: extra frcmod files.
- `gaff`: `0`, `1`, or `2` for no GAFF, GAFF, or GAFF2.
- `group_name`: output prefix, default `MOL`.
- `large_opt`, `lgmodel_chg`, `lgmodel_spin`, `smmodel_chg`, `smmodel_spin`: QM model controls.
- `naa_mol2files`: mol2 files for non-standard amino acids.
- `software_version`: Gaussian/GAMESS flavor.
- `water_model`: `tip3p`, `opc`, etc.

MCPB.py steps:

1. `1`, `1a`, `1m`, `1n`: generate small, standard, and large model PDBs, fingerprint files, residue information, and QM input files.
2. `2`, `2s`, `2e`, `2z`, `2b`: generate pre/final frcmod files. The default `2`/`2s` path uses the Seminario method and needs QM log/fchk output. `2e` is empirical; `2b` creates a blank metal-related frcmod and is useful for smoke tests.
3. `3`, `3a`, `3b`, `3c`, `3d`: fit RESP charges and generate mol2 files for renamed metal-center residues.
4. `4`, `4b`: generate the bonded-model `*_mcpbpy.pdb`, `*_tleap.in`, dry and solvated Amber outputs. `4n1` and `4n2` are nonbonded models and are not the primary target of this task.

For the bonded model, MCPB.py step 4 writes a tleap input that:

- sources the selected protein force field, GAFF/GAFF2, and water model;
- adds atom types for the renamed metal-center residues;
- loads generated `*.mol2` files for the refit residues;
- loads user frcmods and `*_mcpbpy.frcmod`;
- loads `*_mcpbpy.pdb`;
- adds disulfide bonds, metal-ligating bonds, and peptide reconnection bonds for renamed residues;
- saves dry and solvated `prmtop`/`inpcrd`.

## Design Direction

Add an explicit MCPB parameterization path to `AmberParameterizer` instead of encoding MCPB behavior as raw user `additional_tleap_lines`.

The core implementation should introduce a small internal representation for one MCPB run, for example:

```python
MCPBMetalCenterConfig(
    metal: MetalUnit,
    ion_charge: int,
    ion_mol2_path: str | None,
    add_bonded_pairs: list[tuple[Atom, Atom]] | None,
    additional_residues: list[Residue] | None,
    cutoff: float = 2.8,
    force_field: str = "ff19SB",
    gaff: int = 1,
    water_model: str = "opc",
    group_name: str | None = None,
    small_model_charge: int | None = None,
    small_model_spin: int | None = None,
    large_model_charge: int | None = None,
    large_model_spin: int | None = None,
    step2_method: str = "2s",
    step3_method: str = "3b",
    qm_log_path: str | None = None,
    qm_fchk_path: str | None = None,
    premade_outputs_dir: str | None = None,
)
```

The public API can be smaller than this, but these are the facts the implementation must eventually know.

The implementation should keep MCPB.py artifacts in a deterministic subfolder under the parameterizer temp dir, e.g.:

```text
{parameterizer_temp_dir}/mcpb/{group_name}/
```

Expected output object from `_parameterize_metalcenter()`:

```python
MCPBMetalCenterParameter(
    mcpb_pdb_path=".../*_mcpbpy.pdb",
    final_frcmod_path=".../*_mcpbpy.frcmod",
    tleap_in_path=".../*_tleap.in",
    mol2_paths=[".../*.mol2"],
    generated_bond_lines=[...],
    generated_add_atom_types=[...],
    renamed_residue_names=[...],
)
```

The first implementation can parse `*_tleap.in` and reuse its generated `addAtomTypes`, `loadmol2`, `loadamberparams`, and `bond` lines instead of reimplementing MCPB.py's naming logic.

## Development Tasks

### 1. Add MCPB.py Interface Wrapper

Create a small wrapper around MCPB.py, likely in `enzy_htp/_interface/amber_interface.py` initially unless a separate module becomes clearer.

Responsibilities:

- Write MCPB input files from structured options.
- Run MCPB.py steps with clear logging.
- Support at least steps `1`, `2`, `3`, and `4`.
- Allow `step2_method` to be `2b` for cheap smoke tests and `2s` for real Seminario/QM tests.
- Pass `--logf` and `--fchk` for step 2 when required.
- Validate required outputs after each step.

Done when:

- A unit test can write an MCPB input file and verify that all expected lines are present.
- A smoke test can run MCPB.py step `1` on a small fixture and validate that `*_small.pdb`, `*_standard.pdb`, `*_large.pdb`, and fingerprint files exist.

### 2. Generate MCPB Inputs From `Structure`

Teach EnzyHTP how to map a `MetalUnit` and its donor environment to MCPB.py input fields.

Responsibilities:

- Save an MCPB-compatible original PDB with stable atom numbering.
- Resolve `ion_ids` from the saved PDB atom indices, not from stale in-memory indices if renumbering occurred.
- Create or accept `ion_mol2files` for metal ions.
- Convert EnzyHTP `Residue`/`Atom` references to MCPB atom IDs for `add_bonded_pairs`.
- Convert selected extra residues to MCPB `additional_resids`.
- Decide how user charge/spin information is supplied for small and large models.

Done when:

- A test can build an MCPB input file for an existing metal fixture and confirm that `ion_ids`, `original_pdb`, `ion_mol2files`, optional bonded pairs, and optional additional residues match the saved PDB.

### 3. Implement `_parameterize_metalcenter()`

Replace the current TODO in `AmberParameterizer._parameterize_metalcenter()`.

Responsibilities:

- For every configured metal center, run or reuse an MCPB.py output directory.
- Return structured MCPB parameter data instead of only warning.
- Support at least bonded model step `4`/`4b`.
- Keep a cache/reuse mode so expensive QM-derived MCPB outputs can be reused.
- Fail with a clear error when QM-required files for step `2s` are missing.

Done when:

- `AmberParameterizer.run()` no longer warns that metal-center parameterization is unsupported when an MCPB config is supplied.
- Missing MCPB config still gives a clear actionable error or warning, not an obscure tleap failure.

### 4. Integrate MCPB Outputs Into Final tleap

Update `_write_combining_tleap_input()` to consume generated MCPB parameters.

Two viable approaches:

1. Preferred early approach: parse MCPB's generated `*_tleap.in` and transplant only parameter-loading and bond-definition lines into EnzyHTP's normal tleap flow.
2. Later approach: let MCPB.py's `*_mcpbpy.pdb` become the PDB that EnzyHTP loads, then add EnzyHTP solvation/job settings around it.

The implementation must handle:

- `addAtomTypes` blocks.
- `loadmol2` lines for renamed residues.
- `loadamberparams *_mcpbpy.frcmod` and any extra frcmods.
- `bond` lines for metal-donor bonds and peptide reconnection around renamed residues.
- The fact that MCPB.py uses `mol` as the tleap unit name while EnzyHTP currently uses `a`.

Done when:

- The final EnzyHTP-generated tleap input contains the same chemically relevant lines as MCPB.py's manual `*_tleap.in`, normalized for unit name.
- The final `prmtop`/`inpcrd` can be generated by EnzyHTP without manually copying MCPB.py lines into `additional_tleap_lines`.

### 5. Add Workflow Script For End-to-End Validation

Create a new workflow main script, for example:

```text
enzy_htp/workflow_app/mcpb_parameterization_main.py
```

Purpose:

- Demonstrate the intended user-facing MCPB automation.
- Take a prepared metal-containing PDB plus MCPB config/QM artifacts.
- Run EnzyHTP's Amber parameterizer with MCPB enabled.
- Produce final `prmtop`/`inpcrd` and keep comparable intermediate MCPB files.

The script should not be a production workflow template with hard-coded lab paths only. It should be a small reproducible driver that can be pointed at a fixture directory.

Done when:

- Running the script in `new_enzy_htp` produces:
  - MCPB input file;
  - `*_mcpbpy.pdb`;
  - `*_mcpbpy.frcmod`;
  - `*_tleap.in`;
  - EnzyHTP-generated `amber_parm.prmtop`;
  - EnzyHTP-generated `amber_parm.inpcrd`.

### 6. Add Manual-Gold Comparison Fixture

Create a fixture directory with manually generated MCPB.py outputs.

Suggested location:

```text
test/_interface/data/mcpb_manual/
```

It should contain:

- input PDB;
- MCPB input file;
- metal ion mol2 file(s);
- any required QM log/fchk or a documented `2b` smoke-test mode;
- manual `*_mcpbpy.pdb`;
- manual `*_mcpbpy.frcmod`;
- manual `*_tleap.in`;
- manual dry `prmtop`/`inpcrd` if feasible.

Done when:

- A regression test compares EnzyHTP generated MCPB input/tleap content against the manual files after normalizing absolute paths, temporary prefixes, and tleap unit name.
- A higher-cost integration test can optionally compare generated Amber topology against manual topology by checking atom count, residue count, metal-donor bonded pairs, and presence of MCPB atom types/parameters.

## Completion Criteria

This task is complete when all of the following are true.

1. `AmberParameterizer` has an MCPB bonded-model path for metal centers.
2. The MCPB path can be configured without hand-writing final tleap bond lines in `additional_tleap_lines`.
3. Existing ligand and modified-residue parameterization tests still pass.
4. A new workflow script under `enzy_htp/workflow_app/` runs successfully in:

   ```bash
   source ~/bin/miniconda3/bin/activate new_enzy_htp
   ```

5. The workflow script produces EnzyHTP Amber parameters from MCPB.py outputs.
6. The generated MCPB/tleap artifacts are equivalent to a manual MCPB.py run for the fixture case:
   - same MCPB model PDB after normalizing path/order-only differences;
   - same final MCPB frcmod content;
   - same metal-center mol2 residue names and atom types;
   - same metal-donor and peptide reconnection `bond` commands;
   - generated `prmtop`/`inpcrd` pass basic topology checks and contain the metal-center bonded model.
7. The non-MCPB manual path remains available: users can still provide `additional_tleap_lines` directly for special cases.

## Explicit Non-Goals For This Milestone

- Reworking PDB chain-splitting or metal residue categorization.
- Redesigning `MetalUnit` broadly.
- Making protonation or mutation fully metal-aware beyond what is required for MCPB.py input stability.
- Automatically running Gaussian/GAMESS calculations on HPC. The first MCPB implementation may require precomputed QM log/fchk files or use `2b` only for smoke tests.
- Supporting every MCPB.py mode. Bonded model `4`/`4b` is the primary target.

## Risks And Decisions To Make

- Whether to parse MCPB's generated `*_tleap.in` or reimplement its naming and bonding logic. Parsing is lower-risk for the first implementation.
- How to represent MCPB config in the public API without exposing every MCPB.py keyword at once.
- How to cache and reuse expensive QM-derived files cleanly.
- Whether final EnzyHTP parameterization should use MCPB's `*_mcpbpy.pdb` directly or patch EnzyHTP's original saved PDB with MCPB renamed residues.
- How strict topology equivalence should be in CI, given path/order/noise differences in AmberTools output.
