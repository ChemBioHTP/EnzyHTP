import os
from subprocess import SubprocessError
import pytest
import helper

DATA_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/data_dir/"

def test_run_cmd_fail_case():
    cmd_fail = 'cat abc'
    with pytest.raises(SubprocessError) as e:
        helper.run_cmd(cmd_fail, try_time=3, wait_time=1)
    print(e.value)

@pytest.mark.accre
def test_run_cmd_no_retry():
    cmd = 'squeue -u $USER'
    assert len(helper.run_cmd(cmd).stdout) != 0

def test_check_complete_metric_run():
    test_data_path = f"{DATA_DIR}Mutation.dat"
    test_mutants = [
        ['EA323R', 'EB773R', 'GA171R', 'GB621R'],
        ['EA323R', 'EB773R', 'GA171H', 'GB621H'],
        ['EA323R', 'EB773R', 'GA171K', 'GB621K'],
        ['EA323R', 'EB773R', 'GA171D', 'GB621D']]
    result = [helper.check_complete_metric_run(mutant, test_data_path) for mutant in test_mutants]
    assert result == [True, True, False, False]
