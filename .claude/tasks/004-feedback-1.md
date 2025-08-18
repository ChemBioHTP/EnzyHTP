# Feedback for Claude Code

1. For run() and make_job() in alphafold_interface, actually implement it. We need to reference analysis/binding.py. run() or make_job() should only handle the execution of alphafold2 and there should be a outside wrapper predict() to parse the output into EnzyHTP data structures and support the user choice of 1) running alphafold locally (run()) or 2) submit it as a HPC job as a ClusterJob. (make_job())
In detail, in this test environment, alphafold is installed as the colabfold version in the apptainer container. In the command line, it is run by 
```
apptainer run --nv -B ~/bin/colabfold/cache:/cache -B ~/colabfold-test/:/work ~/bin/colabfold_1.5.5-cuda12.2.2.sif   colabfold_batch /work/test.fasta /work/test_out
```
See all options colabfold_batch supports in `./colabfold_options.txt`
2. Make sure you write the code in a way that it supports different Alphafold install. Alphafold mainly have two types of installs, in a docker-like container, or as a python script in the current python environment. The docker-like also can have different types both in the dock-like software (e.g.: docker, apptainer, etc.) and in the alphafold executable and options. (e.g.: the native alphafold using local sequence database or colabfold using mmseq2 server for MSA) Ideally, through class variables in alphafold_config.py. (and may be system_config to support different "docker-like software") An example of alternative install of alphafold runs by:
```
# Set your input/output data path
CALCDIR=/path/to/your/input/and/output/data
# Your input fasta should be in the directory above:
FASTA=CTD-EF.fasta
# Where is the AF2 Inference data
AF2_DATADIR=/sb/apps/alphafold-data.230
# Where is the AF2 Git?
AF2_REPO=/sb/apps/alphafold232/alphafold

python $AF2_REPO/run_alphafold.py \
        --fasta_paths=$FASTA \
        --max_template_date=9999-12-31 \
        --data_dir=$AF2_DATADIR \
        --output_dir=$CALCDIR \
        --uniref90_database_path=$AF2_DATADIR/uniref90/uniref90.fasta \
        --mgnify_database_path=$AF2_DATADIR/mgnify/mgy_clusters_2022_05.fa \
        --uniref30_database_path=$AF2_DATADIR/uniref30/UniRef30_2021_03 \
        --bfd_database_path=$AF2_DATADIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
        --template_mmcif_dir=$AF2_DATADIR/pdb_mmcif/mmcif_files \
        --pdb_seqres_database_path=$AF2_DATADIR/pdb_seqres/pdb_seqres.txt \
        --obsolete_pdbs_path=$AF2_DATADIR/pdb_mmcif/obsolete.dat \
        --uniprot_database_path=$AF2_DATADIR/uniprot/uniprot.fasta \
        --use_gpu_relax
```
Full options about this install is in `./alphafold_run_options.txt`. This way you can get an idea on what may change.
3. Complete unit tests for alphafold_interface.py and prediction.py
