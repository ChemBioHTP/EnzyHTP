# Plan for completing the `residue_pka` module

Please help me complete the development of a residue pKa calculation function in the @enzy_htp/analysis/residue_pka.py. The code there now is very preliminary.

1.  Let's read about other files under @enzy_htp/analysis/, learn about the convention and architecture in this module.
2.  **Create a new `propka_interface.py` module** in `enzy_htp/_interface/` to encapsulate all interactions with the `propka` library. This will include:
    *   A function to run `propka` on a PDB file and return the pKa values.
    *   Error handling for `propka` execution.

3.  **Refactor `enzy_htp/analysis/residue_pka.py`** to:
    *   Define a main "Science API" function `residue_pka()` that takes a `Structure` or `StructureEnsemble` object as input.
    *   Implement a `method` argument to allow for future expansion with other pKa calculation methods (e.g., in the future, we may provide multiple engines such as machine learning based method for users to choose from).
    *   Use a dictionary to map method names to interface functions, similar to the other modules in `enzy_htp/analysis`.
    *   The `residue_pka` function will handle `Structure` and `StructureEnsemble` inputs along with a `Residue` input to specify which residue to calculate pKa for, dispatching the calculation to the appropriate interface function.
    *   Add comprehensive docstrings and type hints to the new functions.
    *   Ensure that the module-level docstring is updated to reflect the new functionality.

4.  **Update `enzy_htp/analysis/__init__.py`** to import and expose the new `residue_pka` function.

5.  **Create a new test file `test/analysis/test_residue_pka.py`** to:
    *   Add unit tests for the `residue_pka` function.
    *   Include tests for both `Structure` and `StructureEnsemble` inputs.
    *   Test the `propka` method.
    *   Use a known protein structure and `propka` results for comparison. (You can leave the test in the state that the function works well but no assert has been made. I will then manually check the correct pKa result and put the answer there.)