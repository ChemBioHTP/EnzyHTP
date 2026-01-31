"""Usage: python af2_main.py <cpu_model_name> <output_work_dir> <num_of_repeats>"""
import glob
import os
import sys, shutil
from enzy_htp.structure_prediction import predict_structure
from enzy_htp.core.clusters.accre_r9 import AccreR9
from enzy_htp.core.job_manager import ClusterJobConfig
from enzy_htp import config as eh_config

def generate_multi_conf_result():
    """Generate multiple conformation results by running multiple predictions on the same sequence with the same MSA input."""
    input_seq = ["VNPTVFFDIAVDGEPLGRVSFELFADKVPKTAENFRALSTGEKGFGYKGSCFHRIIPGFMCQGGDFTRHNGTGGKSIYGEKFEDENFILKHTGPGILSMANAGPNTNGSQFFICTAKTEWLDGKHVVFGKVKEGMNIVEAMERFGSRNGKTSKKITIADCGQL"]
    source_msa_dir = "./af2_output/seq_0/msas"
    num_of_repeat = int(sys.argv[3])  # e.g., 70
    input_seq = input_seq * num_of_repeat # 70 same seq for 70 repeats

    cluster_job_config = ClusterJobConfig(
        cluster = AccreR9(),
        res_keywords = {
            "account": "csb_gpu_acc",
            "partition": "batch_gpu",
            "node_cores": f"{sys.argv[1]}:1",
            "walltime": "3-00:00:00",
            # "account": "yang_lab_csb_iacc",
            # "partition": "interactive_gpu",
            # "qos": "debug_iacc",
            # "node_cores": "nvidia_rtx_a4000:1",
            # "walltime": "30:00",
            }
        )

    eh_config.alphafold.INSTALL_TYPE="alphafold2_native_python"
    eh_config.alphafold.EXECUTABLE_PATH="/sb/apps/alphafold232/alphafold/run_alphafold.py"
    eh_config.alphafold.DATA_DIR="/sb/apps/alphafold-data.230"
    work_dir = sys.argv[2]  # e.g., "./af2_multiconf_output"

    # prepare msa files
    for i in range(num_of_repeat):
        rep_work_dir = f"{work_dir}/seq_{i}"
        rep_msa_dir = f"{rep_work_dir}/msas"
        shutil.copytree(source_msa_dir, rep_msa_dir, dirs_exist_ok=True)

    predict_structure(
        sequences=input_seq, 
        cluster_job_config=cluster_job_config,
        seq_per_job=1,
        array_size=0,
        num_relax=0,
        use_precomputed_msas=True, # TODO support automatically copy files to work_dir/msas/ 
        # precomputed_msa_path_list = ["xxx.a3m"]
        work_dir=work_dir,
        use_templates=False,
        # max_template_date="2021-02-15",
        )

def collect_multi_conf_result():
    """Collect multiple conformation results into a single directory."""
    num_of_repeat = int(sys.argv[3])  # e.g., 70
    work_dir = sys.argv[2]  # e.g., "./af2_multiconf_output"
    collect_dir = f"{work_dir}/collected_results"
    for i in range(num_of_repeat):
        rep_output_dir = f"{work_dir}/seq_{i}"
        src_file = glob.glob(f"{rep_output_dir}/*.pdb")
        dst_dir = f"{collect_dir}/rep_{i}/"
        os.makedirs(dst_dir, exist_ok=True)
        for f in src_file:
            shutil.copy(f, dst_dir)

def main():
    generate_multi_conf_result()
    # collect_multi_conf_result()

if __name__ == "__main__":
    main()