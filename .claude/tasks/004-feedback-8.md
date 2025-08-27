# Feedback for Claude Code

1. In `af2_predict`, different install type also requires different fasta format.
   1) For colabfold install, it use `:` to specify inter-protein chainbreaks for **modeling complexes** (supports homo- and hetro-oligomers). For example **PI...SK:PI...SK** for a homodimer. And it prefers having multiple prediction targets in the same fasta file.
   2) For native alphafold install, it requires each fasta file to only contain one prediction target. If a FASTA file contains multiple sequences, then it will be folded as a multimer. Paths should be separated by commas. All FASTA paths must have a unique basename as the basename is used to name the output directories for each prediction.
   3) To support the multimer input, the `sequences` argument should be in the support both format of [seq1, seq2, seq3, ...] and [[seq1_chain_A, seq1_chain_B, ...], [seq2_chain_A, seq2_chain_B, ...], ...]
2. The arguments and configs related to `_build_alphafold2_native_python_command` needs a refactor. 
   1) `_build_alphafold2_native_python_command` need to support all arguments I put. I also put the related flag in comments. See alphafold_run_options.txt for what they mean.
   2) See the TODO part in the docstring, all these options needs to be supported through alphafold_config.
   3) An important thing is PDB seqres and UniPot are only for model_preset=multimer, error will show up if they are set in a non-multimer run.
3. For these unit tests, you need to have `test_af2_real_accre_r9_native_python_job_setup` pass to claim successfully support the python version of native af2. Keep the test as is and fix any bug. To make the test finish within reasonable time, you need to use precompute msa. I put the one for the test sequence under test/_interface/data/test_a4e23.a3m (The test may take 3 min to finish and you should look for the slurm-xxx.out for the output of the submitted job. Ideally, if that job fails, egg should have a method to report the error information.)
