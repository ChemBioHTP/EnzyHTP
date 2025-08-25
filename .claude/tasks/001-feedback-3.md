# Feedback for Claude Code

Several problems in your code.
    1. Let's write a unit test try to isolate the deepcopy error happened in our lv 5 test.
    2. get_protein_force_field need a unit test. You need to tests different possible input from user. Run the test and make the function more robust. For example, the function is not case sensitive to the input force_field, however should not, this should be addressed. Make sure the test passes.
    3. In run_prepgen, os is already imported, remove the redundant import.
    4. In run_prepgen, the error handling is not robust. Please keep any exception exposed.
    5. The generated frcmod, frcmod2, mol2 file did not use the naming convention. Reference the convention in _parameterize_ligand. (i.e. you need to add target_method in the name.)
    6. We need to make sure the "# 0. search parm lib - same as ligand" part of `_parameterize_modified_res` works. Read and understand how this works and write a seperate unit test to test this. (tip: you need to make the frcmod, frcmod2, mol2 files with the correct name and a new ncaa_lib_... dir under the associated data directory. dont forget the ones for RLP as well)
