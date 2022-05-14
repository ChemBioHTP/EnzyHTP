from subprocess import SubprocessError
import pytest
import helper

def test_run_cmd_fail_case():
    cmd_fail = 'cat abc'
    with pytest.raises(SubprocessError) as e:
        helper.run_cmd(cmd_fail, try_time=3, wait_time=1)

def test_run_cmd_no_retry():
    cmd = 'squeue -u $USER'
    assert len(helper.run_cmd(cmd).stdout) != 0