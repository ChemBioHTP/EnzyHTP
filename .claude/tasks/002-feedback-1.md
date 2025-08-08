# Feedback to Claude Code

1. Do not use nested function in `residue_pka`
2. In the unit test of `residue_pka`, there should not be a if __name__ == __main__ block. The test script should never be called without pytest
3. Run all the unit test. If they fail, fix the problem. (unless it is about the final assert that you would need my manual help) I can see now `test_residue_pka_structure` will not work because the written intermediate PDB have too many solvent residues that it overflow its position in the line. To fix this, we need an option in the PDBParser().get_file_str to not write chain id in the str. If the option is enabled, user also need to receive a warning.
4. The other thing is propka dont really benefit from having solvent in the structure. So after you make sure the test with the solvated structure works. You should add in the get_residue_pka_from_stru to remove solvent by default.