"""Testing enzy_htp.core.cluster.accre.py
Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-09-25
"""
from enzy_htp.core.clusters.accre import Accre

def test_parser_resource_str_gpu():
    res_dict = {
        'core_type' : 'gpu',
        'nodes':'1',
        'node_cores' : '1',
        'job_name' : 'EnzyHTP_MD',
        'partition' : 'pascal',
        'mem_per_core' : '32G',
        'walltime' : '3-00:00:00',
        'account' : 'xxx'
    }
    res_str = Accre().parser_resource_str(res_dict)
    assert res_str == '''#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --job-name=EnzyHTP_MD
#SBATCH --partition=pascal
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --account=xxx
#SBATCH --export=NONE
'''
    res_dict['node_cores'] = '2'
    res_str = Accre().parser_resource_str(res_dict)
    assert res_str == '''#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gres=gpu:2
#SBATCH --job-name=EnzyHTP_MD
#SBATCH --partition=pascal
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --account=xxx
#SBATCH --export=NONE
'''
    res_dict['mem_per_core'] = '10.3G' # round by 0.1
    res_str = Accre().parser_resource_str(res_dict)
    assert res_str == '''#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gres=gpu:2
#SBATCH --job-name=EnzyHTP_MD
#SBATCH --partition=pascal
#SBATCH --mem=21G
#SBATCH --time=3-00:00:00
#SBATCH --account=xxx
#SBATCH --export=NONE
'''