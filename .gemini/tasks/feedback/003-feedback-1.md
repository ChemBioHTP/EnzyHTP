# Feedback to Gemini

1. The selection syntax AmberMask couples with Amber which is not the modular design we want.
2. It also makes no sense to allow user to select freely using the full AmberMask. (e.g.: selecting atoms make that are not full residues make no sense) Let's change the input from a string to a list of residue keys. Dispatch in this API if there are only two element in this list, we treat them as the start and end residue of the domain. (this is a more common use case)
3. You dont need to refactor the original workflow script.
