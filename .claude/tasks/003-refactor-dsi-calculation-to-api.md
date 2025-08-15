# Task: Refactor DSI Calculation Script into an Analysis API

## 1. Analysis of the Existing Script

The current script `enzy_htp/workflow_app/dsi_calculation.py` is a standalone MVP. It correctly calculates the Domain-Domain Interaction Index (DSI) by:
- Taking file paths (`traj_file`, `top_file`) and residue ranges as direct inputs.
- Generating and running a `cpptraj` script to calculate inter-domain distance and individual radii of gyration.
- Parsing the `cpptraj` output to compute the final DSI value for each frame.
- Cleaning up temporary files.

While functional, it's not integrated into the EnzyHTP framework, as it operates on file paths rather than the core `Structure` or `StructureEnsemble` objects.

## 2. Goal
Create a DSI calculation API in the `analysis` module, properly separated from the `amber` interface, that uses standard EnzyHTP data types for selections.

## 3. Plan

### Part 1: Amber Interface Enhancement
A new function will be added to the Amber interface to handle the specific `cpptraj` calculations required for DSI.

- **File:** `enzy_htp/_interface/amber_interface.py`
- **Function name:** `calculate_dsi_metrics`
- **Logic:**
    - This function will encapsulate all `cpptraj` interactions.
    - It also take StructureEnsemble, List[Tuple[str, int]], and List[Tuple[str, int]] as input. (see the dsi api below)
    - It will generate and execute a `cpptraj` script to compute the distance between the centers of mass and the radius of gyration for each of the two domain selections.
    - It will parse the `cpptraj` output, manage all temporary files, and return the dsi result.

### Part 2: New Analysis API
The user-facing API will reside in the `analysis` module and will be decoupled from the Amber implementation details.

- **File:** `enzy_htp/analysis/dsi.py`
- **Function Signature:**
  ```python
  from typing import List, Tuple
  import numpy as np
  from enzy_htp.structure import StructureEnsemble

  def dsi(ensemble: StructureEnsemble, domain1_residues: List[Tuple[str, int]], domain2_residues: List[Tuple[str, int]], engine: str="cpptraj") -> np.ndarray:
      """
      Calculates the Domain-Domain Interaction Index (DSI) for a trajectory.

      DSI is a measure of distance between two domains, defined as:
      DSI = d(com1, com2) - (Rg1 + Rg2)
      where d(com1, com2) is the distance between the centers of mass of the
      two domains, and Rg1 and Rg2 are their respective radii of gyration.

      Args:
          ensemble: A StructureEnsemble object containing topology and trajectory.
          domain1_residues: A list of residue keys (chain_id, residue_idx)
                            for the first domain. If two keys for the same
                            chain are provided, they are treated as the
                            start and end of a continuous residue range.
          domain2_residues: A list of residue keys for the second domain.
          engine: the engine used for the dsi calculation.
      """
  ```
- **Internal Logic:**
    It will use the right interface level function depending on user's choic of engine. (referece other APIs under analysis/ like binding.py)

### Part 3: Development and Testing
1.  **Implement Amber Interface Function:** Add `calculate_dsi_metrics` to `enzy_htp/_interface/amber_interface.py`.
2.  **Implement Analysis API:** Create `enzy_htp/analysis/dsi.py` and implement the `dsi` function
3.  **Add Unit Tests:**
    - Add integration tests for `calculate_dsi_metrics` in `test/_interface/`.
    - In `test/analysis/test_dsi.py`, add unit tests for the `dsi` function. 
