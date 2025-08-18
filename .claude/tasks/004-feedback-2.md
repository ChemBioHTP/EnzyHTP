# Feedback for Claude Code

1. In alphafold_config.py, CONTAINER_BIND_PATHS should not have "/work" or WORK_DIR, it should be part of the logic in the function (run and make_job) in alphafold_interface to empower the `work_dir` argument.
2. Add docstring for each class variables in alphafold_config.py, recycle those tailing comments.
3. In alphafold_interface, for `predict()`:
   1) I renamed it to `af2_predict()` and removed the old `af2_predict()`. This is because there are many versions of alphafold. In your config also make sure to note in varible names those are install for alphafold2 (instead of 1 or 3)
   2) You need to refactor the function for creating fasta file from sequence to a generic helper function in a new file called `sequence.py` in enzy_htp/chemical. Make unit test for it as well. 
   3) You need to put `out_dir` (and all run()/make_job() arguments) explicitly in arguments.
   4) You need to use `fs.clean_temp_file_n_dir(temp_path_list: List[str])` for cleaning files
   5) For the run on cluster case, you need to allow user to specify an `array_size`, `job_check_period` and `seq_per_job` and run it as an array. Instead of `_find_output_files`, you should have make_job() return a `ResultEgg` as a child class of `ModelingResultEgg` to record information of expected output files (see `GaussianQMResultEgg` as an exmaple). Also see `wait_to_array_end_plus` in `enzy_htp/core/job_manager.py` for the api you should call to submit and wait. You are not using ClusterJob correctly now.
4. `subprocess.run` is NOT the function to run a shell command. Use `self.env_manager_.run_command`.
5. You need to have at least one un-mocked unit test for each function you created. Currently may functions you created don't work. Use these un-mocked unit test to fix them!
