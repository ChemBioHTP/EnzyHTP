"""Here is everything job manager need to know about `paracloud-NC-N40` of BSCC (Beijing Super Computer Center.
TODO make a base class for slurm

NOTE: for SLURM based HPC, there will be a common problem with job manager that it creates large
workload to the controller of slurm. You can use scancel -s SIGSTOP to pause the main script to relieve
the workload and use scancel -s SIGCONT to resume the script.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2024-05-20
"""
import re
import os
from subprocess import CompletedProcess, run
import time

from ._interface import ClusterInterface
from ..math_helper import round_by
from ..env_manager import EnvironmentManager
from ..logger import _LOGGER

# Empty env manager just for running commands 
# (env checking will be done whenever a command is needed, otherwise there will be too much warning)
ENV_MANAGER = EnvironmentManager(
    executables=[
        "squeue",
        "sbatch",
        "sacct",
        "scontrol",
        "scancel",
    ]
) # TODO(qz): for some reason ABC is not happy with this going into the class
ENV_MANAGER.check_executables()

class Paratera1(ClusterInterface):
    """
    The CNGrid12 interface
    """
    #############################
    ### External use constant ###
    #############################
    NAME = "CNGrid12"

    # environment presets #
    AMBER_ENV = {
        "CPU": "source /project/wangxin4/soft/amber20_src/amber.sh",
        "GPU": "source /project/wangxin4/soft/amber20_src/amber.sh",
    }

    G16_ENV = {
        "CPU":"""alias g16='/project/wangxin4/soft/gs/g16/g16'
export g16root=/project/wangxin4/soft/gs
export GAUSS_SCRDIR=/project/wangxin4/soft/gs/g16/tmp
source $g16root/g16/bsd/g16.profile"""
    }

    MULTIWFN_ENV = {
        "GPU" : """ulimit -s unlimited
export OMP_STACKSIZE=200M
export Multiwfnpath=/project/wangxin4/soft/Multiwfn_3.8_dev_bin_Linux_noGUI
export PATH=$PATH:/project/wangxin4/soft/Multiwfn_3.8_dev_bin_Linux_noGUI"""
    }

    #############################
    ### Internal use constant ###
    #############################
    # command for submission
    SUBMIT_CMD = "sbatch"
    # regex pattern for extracting job id from stdout
    JOB_ID_PATTERN = r"Submitted batch job ([0-9]+)"
    # command for control & monitor
    KILL_CMD = "scancel"
    HOLD_CMD = "scontrol hold"
    RELEASE_CMD = "scontrol release"
    INFO_CMD = ["squeue", "sacct"] # will check by order if previous one has no info
    # dict of job state
    JOB_STATE_MAP = {
        "pend" : ["CONFIGURING", "PENDING", "REQUEUE_FED", "REQUEUE_HOLD", "REQUEUED"],
        "run" : ["COMPLETING", "RUNNING", "STAGE_OUT"],
        "cancel" : ["CANCELLED", "DEADLINE", "TIMEOUT"],
        "complete" : ["COMPLETED"],
        "error" : ["BOOT_FAIL", "FAILED", "NODE_FAIL", "OUT_OF_MEMORY", "PREEMPTED", "REVOKED", "STOPPED", "SUSPENDED"],
        "exception": ["RESIZING", "SIGNALING", "SPECIAL_EXIT", "RESV_DEL_HOLD" ]
    }

    RES_KEYWORDS_MAP = { 
        "core_type" : None,
        "nodes":"nodes=",
        "node_cores" : {"cpu": "tasks-per-node=", "gpu": "gres=gpu:"},
        "job_name" : "job-name=",
        "partition" : "partition=",
        "mem_per_core" : {"cpu": "mem-per-cpu=", "gpu": "mem="}, # previously using mem-per-gpu= change to mem= (calculate the total memory) base on issue #57
        "walltime" : "time=",
        "clusters" : "clusters=", # TODO this is an example that we may need more conceptual keywords
    }

    PARTITION_VALUE_MAP = {
        "production" : {"cpu": "??",
                        "gpu": "gpu"},
        "debug": {"cpu": "??",
                  "gpu": "gpu"}
    }

    ##########################
    ### Submission Related ###
    ##########################
    @classmethod
    def parser_resource_str(cls, res_dict: dict) -> str:
        """
        1. parser general resource keywords to CNGrid12 specified keywords
        2. format the head of the submission script
        res_dict: the dictionary with general keywords and value
           (Available keys & value format:
                "core_type" : "gpu",
		        "nodes":"1",
                "node_cores" : "5",
                "job_name" : "job_name",
                "partition" : "priv_para",
                "clusters" : "priv",
                "mem_per_core" : "4G",
                "walltime" : "24:00:00",)
        return the string of the resource section
        """
        parsered_res_dict = cls._parser_res_dict(res_dict)
        res_str = cls._format_res_str(parsered_res_dict)
        return res_str

    @classmethod
    def _parser_res_dict(cls, res_dict):
        """
        parser keywords into CNGrid12 style
        """
        new_dict = {}
        for k, v in res_dict.items():
            # remove core type line
            if k == "core_type":
                continue
            # select core type
            if k in ("node_cores", "mem_per_core"):
                new_k = cls.RES_KEYWORDS_MAP[k][res_dict["core_type"]]
            else:
                new_k = cls.RES_KEYWORDS_MAP[k]
            new_dict[new_k] = v
        # process total mem
        for k in new_dict:
            if k == "mem=":
                mem_per_core_n_gb = new_dict[k].rstrip("GB")
                total_mem = round_by(float(mem_per_core_n_gb) * float(res_dict["node_cores"]), 0.1) # round up
                new_dict[k] = f"{total_mem}G"
        return new_dict

    @staticmethod
    def _format_res_str(parsered_res_dict):
        """
        format parsered dictionary
        NOTE export=NONE is necessary to keep env decoupled
        """
        res_str = "#!/bin/bash\n"
        for k, v in parsered_res_dict.items():
            res_line = f"#SBATCH --{k}{v}\n"
            res_str += res_line
        # res_str += "#SBATCH --clusters=swarm"
        return res_str
    
    @classmethod
    def submit_job(cls, sub_dir, script_path, debug=0) -> tuple[str, str]:
        """
        submit job submission script to the cluster queue
        cd to the sub_path and run submit cmd
        Return:
            (job_id, slurm_log_file_path)
        Raise:
            SubprocessError if sbatch wont work after 1440 retries (12hrs)

        CNGrid12 sbatch rule:
            stdout: Submitted batch job ######## on cluster XXX
            file: slurm-#######.out will be generated in the *submission dir*
            exec: commands in the script will run under the *submission dir*
        """
        # for unit test
        if debug:
            return (f"{cls.SUBMIT_CMD} {script_path}", sub_dir, script_path), None

        cwd = os.getcwd()
        # cd to sub_path
        os.chdir(sub_dir)
        try:
            submit_cmd = ENV_MANAGER.run_command(exe=cls.SUBMIT_CMD,
                                                 args=[script_path],
                                                 try_time=1440, wait_time=60, timeout=120,
                                                 log_level="debug") # 12 hrs
        finally:
            # TODO(shaoqz) timeout condition is hard to test
            os.chdir(cwd) # avoid messing up the dir
        
        job_id = cls._get_job_id_from_submit(submit_cmd)
        slurm_log_path = cls._get_log_from_id(sub_dir, job_id)
        return (job_id, slurm_log_path)

    @classmethod
    def _get_job_id_from_submit(cls, submit_job: CompletedProcess) -> str:
        """
        extract job id from the output of the submission command run.
        """
        job_id = re.search(cls.JOB_ID_PATTERN, submit_job.stdout).group(1)
        return job_id

    @classmethod
    def _get_log_from_id(cls, sub_dir: str, job_id: str) -> str:
        """
        file: slurm-#######.out will be generated in the *submission dir*
        """
        file_path = sub_dir + f"/slurm-{job_id}.out"
        return file_path

    ###############################
    ### Post-submission Related ###
    ###############################

    @classmethod
    def kill_job(cls, job_id: str) -> CompletedProcess:
        cmd = f"{cls.KILL_CMD} {job_id}"
        kill_cmd = run(cmd, timeout=20, check=True,  text=True, shell=True, capture_output=True)
        return kill_cmd

    @classmethod
    def hold_job(cls, job_id: str) -> CompletedProcess:
        cmd = f"{cls.HOLD_CMD} {job_id}"
        hold_cmd = run(cmd, timeout=20, check=True,  text=True, shell=True, capture_output=True)
        return hold_cmd

    @classmethod
    def release_job(cls, job_id: str) -> CompletedProcess:
        cmd = f"{cls.RELEASE_CMD} {job_id}"
        release_cmd = run(cmd, timeout=20, check=True,  text=True, shell=True, capture_output=True)
        return release_cmd

    @classmethod
    def get_job_info(cls, job_id: str, field: str, wait_time=3) -> str:
        """
        get information about the job_id job by field keyword
        1. use squeue frist (fast) and
        if nothing is found (job finished)
        wait for wait time and 2. use sacct (slow)
        * will try each command (if fail) up to 2880 times with a 30s gap (covers 1 day time span)
        Arg:
            job_id
            field: supported keywords can be found at https://slurm.schedmd.com/sacct.html *can only take one keyword at a time*
            wait_time: for sacct run in second (default: 3s)
            **The `sacct` command takes some time (1-5s) to update the information of the job**
        Return:
            The field value as a string
        """
        # get info
        # squeue
        cmd = f"{cls.INFO_CMD[0]} -u $USER -O JobID,{field}" # donot use the -j method to be more stable
        info_run = ENV_MANAGER.run_command(exe=cmd.split()[0],
                                           args=cmd.split()[1:],
                                           try_time=2880, wait_time=120, timeout=120,
                                           log_level="debug")
        # if exist
        info_out_lines = info_run.stdout.strip().splitlines()
        for info_line in info_out_lines: 
            info_line_parts = info_line.strip().split()
            if len(info_line_parts) < 2:
                _LOGGER.debug(f"field: {field} is not supported in squeue. Switch to sacct.")
                break
            if job_id in info_line:
                job_field_info = info_line_parts[1].strip().strip("+")
                return job_field_info
        # use sacct if squeue do not have info
        # wait a update gap
        _LOGGER.debug("No info from squeue. Switch to sacct")
        time.sleep(wait_time)
        cmd = f"{cls.INFO_CMD[1]} -j {job_id} -o {field}"
        info_run = ENV_MANAGER.run_command(exe=cmd.split()[0],
                                           args=cmd.split()[1:],
                                           try_time=2880, wait_time=120, timeout=120,
                                           log_level="debug")          
        # if exist
        info_out = info_run.stdout.strip().splitlines()
        if len(info_out) >= 3:
            job_field_info = info_out[2].strip().strip("+")
            return job_field_info
        raise Exception(f"No information is found for {job_id}")
    
    @classmethod
    def get_job_state(cls, job_id: str) -> tuple[str, str]:
        """
        determine if the job is:
        Pend or Run or Complete or Canel or Error
        Return: 
            a tuple of
            (a str of pend or run or complete or canel or error,
                the real keyword form the cluster)
        """
        state = cls.get_job_info(job_id, "State")
        for k, v in cls.JOB_STATE_MAP.items():
            if state in v:
                return (k, state)
        raise Exception(f"Do not regonize state: {state}")
