# Feedback for Gemini

0. I made some changes to the code you wrote. Be aware of it.
1. The prediction.py need a dictionary to contain all supported engines. See analysis/dsi.py for example.
2. For run() in alphafold_interface actually implement it. We need to reference analysis/binding.py. run() should only handle the execution of alphafold2 and there should be a outside wrapper predict() to parse the output into EnzyHTP data structures and support the user choice of 1) running alphafold locally or 2) submit it as a HPC job as a ClusterJob.
You can reference the following scratch code for the run() method, fit them in EnzyHTP's design principle:
    ```
    def run(
        sequences: Union[str, List[str]],
        out_dir: Union[str, Path],
        database_source: str,     # path to your database directory or use colabfold and mmseq2 in this case equal to the msa_mode option
        num_models: int,
        num_recycles: int,
        num_relax: int = 0,
        relax_max_iteration: int = 200,
        use_templates: bool = False,
        max_template_date: Optional[str] = None,
        model_preset: Optional[str] = None,   # monomer/monomer_ptm/multimer
        db_preset: str = "reduced_dbs",       # or "full_dbs"
        addition_options: Optional[List[str]] = None,
    ) -> Dict[str, str]:
    # the config module controls the actual command of AF2. (either it is a python script or an docker/apptainer command; either use colabfold with mmseq2 server or alphafold with local database)
    ```
3. You also need unit tests for alphafold_interface.
