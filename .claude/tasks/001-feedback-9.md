# Feedback for Claude Code

Several problems in your code.
    1. `run_prepgen` needs a rework. I have already figured out the bug we faced previously about the path length of mc file. I summarized my note in the docstring. The plan is:
       a. add an argument: work_dir, by default it is the current dir, we will change to this dir to run prepgen because it generates many intermediate files under work dir
       b. in the `_parameterize_modified_res` function, we need to specify a work dir, let be `self.parameterizer_temp_dir`.
       c. check the length of the absolute mc path. If it >20 characters, we need to copy it over to the working dir and use a very short name `temp.mc`. Make sure you use `fs.get_valid_temp_name` to avoid overwrite.
       d. the current code also dont clean up the intermediate files correctly. They need to be cleaned before changing the dir back.
       e. finally run unit tests about prepgen and make sure the code works well.
    2. The `test_amber_parameterizer_run_lv_5` still don't work. The generated prmtop file is empty and there is no inpcrd file generated. I commented out the clean up line in this test so that I can check your result. Your work is not done unless the test passes.
    TIPS: **Primary Assertion:** Assert that the final `.prmtop` and `.inpcrd` files are created and are not empty by calling `params.is_valid()`.

!! RUN unit tests related to this problem. Create unit test when it is necessary. A problem is only considered resolved when the related unit test passes.

 
