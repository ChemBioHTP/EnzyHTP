# Feedback for Claude Code

1. For `AlphaFold2ResultEgg`, there should only be one job associated with one result egg. Change your design accordingly.
2. `_find_output_files` should just find return filename map to the file path. There should be another function filter and map sequence to a group of files and screen for the best file out of this group for this seqence (i.e.: highest rand and relaxed)
3. The final return of `af2_predict` should have the full seqence itself as key.
4. In `prediction.py`, sequences should as accept a path of fasta file which is then parsed to a list of seqences in the function.