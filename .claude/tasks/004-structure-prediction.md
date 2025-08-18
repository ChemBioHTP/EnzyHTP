# Task: Add Structure Prediction Science API

This task is to add a new Science API module for predicting protein structure from an amino acid sequence.

## 1. Create New Module and API Files

-   Create a new directory: `enzy_htp/structure_prediction/`
-   Create `enzy_htp/structure_prediction/__init__.py`
-   Create the main API file: `enzy_htp/structure_prediction/prediction.py`

## 2. Design and Implement the Science API

-   In `enzy_htp/structure_prediction/prediction.py`, define a function `predict_structure`.
-   `predict_structure` will take arguments like `sequence`, `engine='alphafold'`, and other relevant parameters.
-   The function will return a `Structure` object.
-   The function will use the interface pattern to call the specified engine.

## 3. Create AlphaFold Interface

-   Create a new config file: `enzy_htp/_config/alphafold_config.py`. This will contain the configuration for AlphaFold (e.g., path to executable, database paths).
-   Create a new interface file: `enzy_htp/_interface/alphafold_interface.py`.
-   The interface will have a class, e.g., `Alphafold`, with a method like `predict(sequence)` that executes the AlphaFold prediction and returns the path to the predicted PDB file.
-   Update `enzy_htp/_config/config.py` to include the new `alphafold_config`.
-   Update `enzy_htp/_interface/interface.py` to include the new `alphafold_interface`.

## 4. Add Tests

-   Create a new test directory: `test/structure_prediction/`
-   Create `test/structure_prediction/test_prediction.py`.
-   Add unit tests for `predict_structure`. This will likely involve mocking the AlphaFold interface to avoid running actual predictions during testing.
-   If possible, add an integration test that runs a small AlphaFold prediction.

## 5. Documentation

-   Add comprehensive docstrings to all new modules, classes, and functions, explaining their purpose, parameters, and return values.
-   Update the project's main `README.md` or other relevant documentation to include information about the new structure prediction module.

## 6. Refinements

-   Ensure the new module is properly integrated into the `enzy_htp` package by updating `__init__.py` files.
-   Make sure the code follows the project's coding style and conventions.
