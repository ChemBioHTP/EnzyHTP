# Feedback for Claude Code

1. `af2_predict` need to have a core type option that reflect on commands, env when use cluster job, etc.
2. We need default values for the res_keywords for AF2 cluster job. Reference how it is done in `get_mmpbgbsa_energy` from `amber_interface.py`.
3. We need unit test in `test/_interface/test_alphafold_interface.py` for the alternative install of alphafold2 and HPC submission. I have put you on a HPC environment (ACCRE R9), where you can test the AF2 install and when cluster_job_config is not None. A scratch submission script for running the AF2 install here is (use this as a reference only, we dont like how this script is designed):
   ```
    #SBATCH --account=csb_gpu_acc
    #SBATCH --partition=batch_gpu
    #SBATCH --constraint=csbtmp
    #SBATCH --mail-type=BEGIN,END,FAIL
    #SBATCH --nodes=1
    #SBATCH --ntasks=6
    #SBATCH --gres=gpu:nvidia_rtx_a6000:1
    #SBATCH --mem=24G
    #SBATCH --time=16:00:00
    #SBATCH --job-name=af232-test
    #SBATCH --output=af232-test.log

    # Set your input/output data path
    CALCDIR=/path/to/your/input/and/output/data

    # Your input fasta should be in the directory above:
    FASTA=CTD-EF.fasta
    # Where is the AF2 miniconda environment
    AF2_MINICONDA=/sb/apps/alphafold232/miniconda3
    # Where is the AF2 Inference data
    AF2_DATADIR=/sb/apps/alphafold-data.230
    # Where is the AF2 Git?
    AF2_REPO=/sb/apps/alphafold232/alphafold

    cd $CALCDIR

    #Look at the driver and GPUs
    nvidia-smi

    echo -n "Running on "
    echo $SLURM_JOB_NODELIST

    # Activate CSB Alphafold2 miniconda environment
    source $AF2_MINICONDA/bin/activate af232
    export LD_LIBRARY_PATH=$AF2_MINICONDA/envs/af232/lib:$LD_LIBRARY_PATH

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
4. For these unit tests:
   1) see `test_amber_md_step_make_job_new_accre_node_cores` for how the cluster_job_config is set up on Accre R9.
   2) In your test, change the values in AlphafoldConfig temp to match in install on Accre R9. The install type should be alphafold native python. See more from the script example in comment 3.
   3) For tests that you will actually submit jobs, use:
   ```
   account : yang_lab_csb_iacc
   partition : interactive_gpu
   qos : debug_iacc
   node_cores : nvidia_rtx_a4000:1
   ```
   4) We need at least one test that does not mock anything.