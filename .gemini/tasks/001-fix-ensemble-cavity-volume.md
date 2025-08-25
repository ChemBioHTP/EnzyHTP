# Task: Fix and Enhance `ensemble_cavity_volumes`

## 1. Objective

The primary goal of this task is to fix the `ensemble_cavity_volumes` function located in `enzy_htp/analysis/cavity.py`, which is currently not functioning as expected. The immediate measure of success is to make the existing unit test, `test_ensemble_cavity_volumes`, pass.

Furthermore, the task includes expanding the test coverage for `ensemble_cavity_volumes` to ensure its robustness and correctness across all its main functionalities.

## 2. Problem Description

The `ensemble_cavity_volumes` function is intended to calculate the volume of a specific cavity across all structures in a `StructureEnsemble`. However, it is currently failing its unit test `test_ensemble_cavity_volumes`.

The function has three distinct modes for identifying the target cavity to be tracked across the ensemble:
1.  By a list of composing residues (`composing_residues`).
2.  By a ligand selection string (`contain_ligand`).
3.  By a reference `Cavity` object (`target_cavity`).

## 3. Scope of Work

### Phase 1: Bug Fix

1.  **Analyze the bug:** Run `test_ensemble_cavity_volumes` first. Analyze the error information. Investigate the `ensemble_cavity_volumes` function in `enzy_htp/analysis/cavity.py` and the corresponding test `test_ensemble_cavity_volumes` to understand the root cause of the failure.
2.  **Implement the fix:** Correct the logic within `ensemble_cavity_volumes` to ensure it correctly identifies and tracks the cavity across the ensemble frames for all specification modes.
3.  **Verify the fix:** Run the `test_ensemble_cavity_volumes` test and ensure it passes.

### Phase 2: Test Coverage Expansion

Once the function is fixed, the next step is to improve its test coverage.

1.  **Review existing tests:** Examine the current `test_ensemble_cavity_volumes` to understand what it covers.
2.  **Implement new unit tests:** Add new tests to specifically and independently validate each of the three cavity selection mechanisms:
    *   A test case for cavity selection using the `composing_residues` argument.
    *   A test case for cavity selection using the `contain_ligand` argument.
    *   A test case for cavity selection using the `target_cavity` argument.

These tests should be comprehensive and check for edge cases where applicable.

## 4. Relevant Files

*   **Function to be fixed:** `/panfs/accrepfs.vampire/home/shaoq1/bin/EnzyHTP-ref/enzy_htp/analysis/cavity.py`
*   **Existing unit test:** The location of `test_ensemble_cavity_volumes` needs to be identified. It is likely in `test/analysis/test_cavity.py`.
*   **File for new unit tests:** Same as the existing unit test file.

## 5. Definition of Done

*   The `ensemble_cavity_volumes` function is fixed and works as intended.
*   The original `test_ensemble_cavity_volumes` test passes.
*   New, separate unit tests for each of the three cavity specification modes (`composing_residues`, `contain_ligand`, `target_cavity`) are implemented and pass.
*   The code adheres to the project's coding style and conventions.
