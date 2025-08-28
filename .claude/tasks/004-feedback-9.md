# Feedback for Claude Code

1. In `af2_predict`, different install type also requires different fasta format.
   1) For colabfold install, it use `:` to specify inter-protein chainbreaks for **modeling complexes** (supports homo- and hetro-oligomers). For example **PI...SK:PI...SK** for a homodimer. And it prefers having multiple prediction targets in the same fasta file.
   2) For native alphafold install, it requires each fasta file to only contain one prediction target. If a FASTA file contains multiple sequences, then it will be folded as a multimer. Paths should be separated by commas. All FASTA paths must have a unique basename as the basename is used to name the output directories for each prediction.
   3) To support the multimer input, the `sequences` argument should be in the support both format of [seq1, seq2, seq3, ...] and [[seq1_chain_A, seq1_chain_B, ...], [seq2_chain_A, seq2_chain_B, ...], ...]
   
   Currently, the support for multimer input when native alphafold install is not working. In native alphafold install's case, a multimer containing input should create multiple fasta files each containing a prediction target. (The currrent code only creates one fasta file)
2. In `_build_alphafold2_native_python_command`, please create and use functions to obtain database paths and exe paths from the alphafold_config module instead of using directly the class variable. That is, you need to create function in alphafold_config for getting each db/exe path, so that unless the corresponding class variable is explicitly defined, it will grow the path from the DATA_DIR.
   For example:
   In alphafold_config.py:
      DATA_DIR = "/sb/apps/alphafold-data.230"
      UNIREF90_DATABASE_PATH = "{DATA_DIR}/uniref90/uniref90.fasta"

      def get_uniref90_database_path(self):
         return UNIREF90_DATABASE_PATH.format(DATA_DIR=self.DATA_DIR)

