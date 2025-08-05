# Task: Implement Modified Amino Acid (MAA) Support in AmberParameterizer

## Objective

The goal of this task is to complete the implementation of the `_parameterize_modified_res` method in the `AmberParameterizer` class, located in `enzy_htp/_interface/amber_interface.py`. This will enable the parameterization of structures containing modified amino acids for Amber MD simulations.

## Background

The `AmberParameterizer` is a key component of EnzyHTP, responsible for generating the necessary topology and coordinate files for running molecular dynamics simulations with Amber. While the class has stubs for handling modified amino acids (MAAs), the implementation is currently incomplete. The existing code raises a `NotImplementedError` and contains placeholder comments outlining the required steps. This task will fill in the missing logic to make MAA parameterization functional.

## Requirements

1.  **Complete the `_parameterize_modified_res` method:** Implement the full workflow for parameterizing a single modified amino acid. This includes:
    *   Generating an `.ac` file using `antechamber`.
    *   Fixing atom types in the `.ac` file.
    *   Creating a `.mc` file.
    *   Running `prepgen` to produce a `.prepin` file.
    *   Using `antechamber` again to get a `.mol2` file.
    *   Running `parmchk2` to generate the final `.frcmod` files.
2.  **Handle file I/O:** Ensure that all intermediate and final parameter files are created in the appropriate locations, using the existing file system utilities (`fs`).
3.  **Return correct values:** The method should return a tuple containing the path to the generated `mol2` file and a list of paths to the `frcmod` files.
4.  **Integration:** Ensure the completed method integrates correctly with the `AmberParameterizer.run` method.

## Implementation Plan

The agent should focus its work on the `_parameterize_modified_res` method within `enzy_htp/_interface/amber_interface.py`.

1.  **Familiarize with existing code:** Review the `_parameterize_ligand` method to understand the existing half-done code for parameterizing non-canonical residues.
2.  **Implement `antechamber` call:** Use `self.parent_interface.antechamber_ncaa_to_moldesc` to generate the initial `.ac` file for the MAA.
3.  **Address atom type correction:** Implement the logic to correct atom types. This may require parsing the `.ac` file and applying corrections based on the MAA's structure.
4.  **Implement `make_mc_file` call:** Use `self.parent_interface.make_mc_file` to generate the `.mc` file.
5.  **Implement `prepgen` call:** A new helper method in `AmberInterface` might be needed to wrap the `prepgen` command-line tool.
6.  **Implement second `antechamber` call:** Use `self.parent_interface.run_antechamber` to convert the `.prepin` file to a `.mol2` file.
7.  **Implement `parmchk2` calls:** Use `self.parent_interface.run_parmchk2` to generate the required `.frcmod` files.
8.  **Clean-up:** Remove the `raise Exception("TODO")` and ensure all temporary files are handled correctly.

## Testing Plan

To verify the implementation, the agent should:

1.  **Create a new test case:** Add a new test function in `test/interface/test_amber_interface.py`.
2.  **Use a real MAA:** The test should use a `Structure` object containing a known modified amino acid (e.g., phosphotyrosine).
3.  **Run the parameterizer:** The test should instantiate `AmberParameterizer` and call the `run` method on the test structure.
4.  **Assert file creation:** The test should assert that the final `.inpcrd` and `.prmtop` files are created and are not empty.
5.  **(Optional) Inspect generated files:** For a more thorough test, the generated `.mol2` and `.frcmod` files for the MAA can be inspected to ensure they are chemically reasonable.

## Definition of Done

- [ ] The `_parameterize_modified_res` method in `enzy_htp/_interface/amber_interface.py` is fully implemented.
- [ ] The method no longer raises a `NotImplementedError` or `Exception`.
- [ ] The method correctly generates `.mol2` and `.frcmod` files for a given MAA.
- [ ] A new unit test that successfully parameterizes a structure with an MAA has been added and passes.
