# Feedback for Claude Code

Several problems in your code.
    1. the unit test for test_amber_parameterizer_run_lv_5 is not passed. You should read the error information and debug until it passes.
    2. In amber_interface.py, the logic for finding the protein force field should be wrapped as a method. A list of supported force field should be stored as a class variable. (reference other mappers in AmberInterface)
    3. ff19SB use XC for CA, I already made edit to the file and fixed it for you.
