# Task: Refactor DSI Calculation Script into an Analysis API

## 1. Analysis of the Existing Script

The current script `enzy_htp/workflow_app/dsi_calculation.py` is a standalone MVP. It correctly calculates the Domain-Domain Interaction Index (DSI) by:
- Taking file paths (`traj_file`, `top_file`) and residue ranges as direct inputs.
- Generating and running a `cpptraj` script to calculate inter-domain distance and individual radii of gyration.
- Parsing the `cpptraj` output to compute the final DSI value for each frame.
- Cleaning up temporary files.

While functional, it's not integrated into the EnzyHTP framework, as it operates on file paths rather than the core `Structure` or `StructureEnsemble` objects.

## 2. Proposed API Design

A new API function will be created within the `analysis` module.

- **New File:** `enzy_htp/analysis/dsi.py`
- **Function Signature:**
  ```python
  import numpy as np
  from enzy_htp.structure import StructureEnsemble

  def dsi(ensemble: StructureEnsemble, selection1: str, selection2: str) -> np.ndarray:
      """
      Calculates the Domain-Domain Interaction Index (DSI) for a trajectory.

      DSI is a measure of the interaction between two domains, defined as:
      DSI = d(com1, com2) - (Rg1 + Rg2)
      where d(com1, com2) is the distance between the centers of mass of the
      two domains, and Rg1 and Rg2 are their respective radii of gyration.

      Args:
          ensemble: A StructureEnsemble object containing the system's
                    topology and trajectory.
          selection1: An AmberMask selection string for the first domain.
                      (e.g., ":1-418&!@H=").
          selection2: An AmberMask selection string for the second domain.
                      (e.g., ":419-515&!@H=").

      Returns:
          A numpy array of DSI values, one for each frame in the trajectory.
      """
  ```
- **Design Rationale:**
    - **Input:** The primary input is a `StructureEnsemble` object. This aligns with EnzyHTP's core design principle of being "Structure-Centric". The `StructureEnsemble` encapsulates the topology and trajectory data, abstracting away file paths from the user. The domains are specified using AmberMask selection strings, which is a flexible and powerful way to define molecular regions and is native to the `cpptraj` interface being used.
    - **Output:** The function will return a `numpy.ndarray` of DSI values. This is a standard and convenient format for any subsequent scientific analysis or plotting.

## 3. Development Plan

The following steps will be taken to implement the new API:

1.  **Create the API Module:**
    - Create the new file: `enzy_htp/analysis/dsi.py`.

2.  **Implement the `dsi` function:**
    - Define the function with the signature proposed above.
    - Adapt the core logic from the existing script. This involves:
        - Getting the necessary topology and trajectory file paths from the input `ensemble` object.
        - Using the `selection1` and `selection2` arguments directly in the `cpptraj` input script.
        - Leveraging `enzy_htp.core.file_system` for robust temporary file management.
        - Parsing the output and returning the results as a NumPy array.

3.  **Refactor the Original Workflow Script:**
    - Modify `enzy_htp/workflow_app/dsi_calculation.py`.
    - Remove the now-redundant `dsi` function from this file.
    - Update its `main` function to use the new API. This will involve creating a `StructureEnsemble` object from the sample data and calling `enzy_htp.analysis.dsi.dsi`, making the script a clean example of how to use the new API.

4.  **Add Unit Tests:**
    - Create a new test file: `test/analysis/test_dsi.py`.
    - Add a test that:
        1. Creates or loads a sample `StructureEnsemble` from the test data.
        2. Calls the `dsi` function with known selections.
        3. Asserts that the output is a `np.ndarray` with the correct shape and that its values are consistent with a pre-calculated result.
