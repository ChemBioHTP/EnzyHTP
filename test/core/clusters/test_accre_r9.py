"""Testing enzy_htp.core.cluster.accre_r9.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2025-03-19
"""
from enzy_htp.core.clusters.accre_r9 import AccreR9

def test_parser_resource_str_gpu():
    res_dict = {
        'core_type' : 'gpu',
        'nodes':'1',
        'node_cores' : 'nvidia_rtx_a6000:1',
        'job_name' : 'EnzyHTP_MD',
        'partition' : 'batch_gpu',
        'mem_per_core' : '32G',
        'walltime' : '3-00:00:00',
        'account' : 'xxx'
    }
    res_str = AccreR9().parser_resource_str(res_dict)
    assert res_str == '''#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gres=gpu:nvidia_rtx_a6000:1
#SBATCH --job-name=EnzyHTP_MD
#SBATCH --partition=batch_gpu
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --account=xxx
#SBATCH --export=NONE
'''
    res_dict['node_cores'] = 'nvidia_rtx_a6000:2'
    res_str = AccreR9().parser_resource_str(res_dict)
    assert res_str == '''#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gres=gpu:nvidia_rtx_a6000:2
#SBATCH --job-name=EnzyHTP_MD
#SBATCH --partition=batch_gpu
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --account=xxx
#SBATCH --export=NONE
'''
    res_dict['mem_per_core'] = '10.3G' # round by 0.1
    res_str = AccreR9().parser_resource_str(res_dict)
    assert res_str == '''#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gres=gpu:nvidia_rtx_a6000:2
#SBATCH --job-name=EnzyHTP_MD
#SBATCH --partition=batch_gpu
#SBATCH --mem=21G
#SBATCH --time=3-00:00:00
#SBATCH --account=xxx
#SBATCH --export=NONE
'''
