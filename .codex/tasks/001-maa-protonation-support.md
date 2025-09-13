# Task: Modified Amino Acid Protonation Support in `protonate_stru`

Owner: preparation/protonation
Status: Planning

## Goals
- Keep mode: allow users to preserve existing protonation of modified amino acids (mod-AAs) while protonating canonical residues via PDB2PQR.
- Auto mode: provide an autonomous path to protonate mod-AAs robustly without breaking peptide connectivity.

## Scope
- Only alter mod-AA protonation behavior in `enzy_htp/preparation/protonate.py::protonate_stru` and related helpers.
- Preserve current flows:
  - Peptide (canonical AAs): PDB2PQR (initial library + PROPKA3).
  - Ligands: PyBel when requested.
- Do not change other modules’ behavior beyond what’s required for mod-AA support.

## API Changes (Backward Compatible)
- Extend `protonate_stru` parameters:
  - `mod_aa_strategy: Literal['keep','pybel'] = 'keep'`
  - Optional `mod_aa_engine` alias (deprecated synonym) for compatibility, if needed.
- Document strategies and caveats in the function docstring.
- Log the selected strategy.

## Behavior Summary
- Default (`mod_aa_strategy='keep'`):
  - Write full structure to PDB2PQR so PROPKA3 can consider ligands; update only canonical residues back into the structure using existing peptide-only update (`remove_non_peptide`, `clone_residue_keys(amino_acid_only=True)`, `update_residues`). Mod-AAs remain unchanged.
- `pybel`:
  - Per mod-AA residue: generate a capped fragment (N-term H, C-term OH); use PyBel at given pH to generate hydrogens; merge sidechain hydrogens back while protecting backbone heavy atoms and their peptide-bond connectivity.

## Implementation Plan

1) Characterize PDB2PQR on mod-AAs
- Add a small test fixture structure containing canonical peptide plus one `ModifiedResidue`.
- Empirically confirm: PDB2PQR output omits non-library mod-AA or leaves them unchanged; current peptide-only update does not alter mod-AA.
- Confirm `clone_residue_keys(amino_acid_only=True)` keeps canonical mapping valid when mod-AAs are ignored.

2) Add `mod_aa_strategy` API
- Update `protonate_stru` signature and docstring.
- Parse and validate strategy; default to `keep`.
- Emit informative logs about selected behavior.

3) KEEP mode
- No changes to canonical peptide pipeline.
- After PDB2PQR, proceed with existing peptide-only update; skip any operation on mod-AAs.
- Test: verify canonical residues gain hydrogens while mod-AA is byte-for-byte identical (atom names, count, coordinates and hydrogens).

4) PyBel-based mod-AA protonation
- For each `ModifiedResidue`:
  - Export residue PDB (without H) using `PDBParser.get_file_str(residue)`.
  - Run `pybel_protonate_pdb_ligand()` at target pH; provide `ref_name_path` to `_fix_pybel_output` to preserve heavy-atom names.
  - Read back the protonated residue.
  - Merge:
    - Identify and protect backbone atoms.
    - Remove original sidechain hydrogens.
    - From PyBel result, collect hydrogens attached to sidechain heavy atoms (name-based mapping), add to original residue.
  - Reindex atoms and maintain parent links; handle errors by logging and skipping the residue.

5) Integrate and log
- Wire mod-AA branch after peptide update and before metal donor fixes.
- Clear, concise logs:
  - Strategy selection, per-residue success/failure, fallbacks taken, and any skipped residues.
- Ensure temp files cleanup mirrors existing ligand flow.

6) Tests (pytest; targeted only)
- KEEP mode: canonical residues change; mod-AA unchanged (atom set, H count, coordinates).
- PyBel mode: simple mod-AA missing sidechain H gains hydrogens; backbone unchanged.
- Antechamber mode: conditional test (skip if AmberTools not available); residue becomes hydrogen-complete and peptide connectivity intact.
- PDB2PQR behavior: confirm mod-AA not updated by peptide-only update; residue-key cloning remains valid.

7) Docs and examples
- Update `protonate_stru` docstring and prep docs:
  - Explain `mod_aa_strategy` options and trade-offs (pH-aware vs robustness).
  - Note antechamber requirements (charge/spin) and non-pH-aware nature.
  - Provide minimal usage examples for KEEP and AUTO modes.

## Risks and Mitigations
- Name alignment failures for complex mod-AAs in PyBel output: mitigate by `_fix_pybel_output` and name-based mapping; otherwise log and skip.
- Connectivity mishaps when replacing residues: rely on existing `init_connectivity` + guardrails; add unit tests around peptide bond integrity.
