# Feedback for Claude Code

1. `make_job` need to keep up with the new multimer logic now when distributing jobs into the array.
2. The current logic of passing `sequences` to the fasta file and final cmd line involves many duplicated code. The sequence is first converted to fasta file and then back to sequence and then into fasta file again. Please refactor related logic in `af2_predict`, `make_job`, `run` so that the sequences are only converted to fasta file(s) once and the multimer related logics is handle is a clean and readable way.
