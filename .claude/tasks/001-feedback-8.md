# Feedback for Claude Code

Several problems in your code.
    1. `test_run_prepgen` is almost there. The only thing to fix is the path from the input arguments are dealt in a messy way right now. Is it possible that mc file and in file are not in the same directory. Using basename is not a good idea. You should use `_reduce_path_in_mdin` in AmberInterface instead. Try it and test it with all prepgen tests. And fix any possible problem.
    2. In `test_run_prepgen_failure_cases`, you should add another case when the ac file is nonsense. Store this file under the data dir of the test folder of _interface and use it in the test.
    3. The `test_amber_parameterizer_run_lv_5` still don't work. The generated prmtop file is empty and there is no inpcrd file generated. I commented out the clean up line in this test so that I can check your result. Work until the test works. 
    TIPS: **Primary Assertion:** Assert that the final `.prmtop` and `.inpcrd` files are created and are not empty by calling `params.is_valid()`.

!! RUN unit tests related to this problem. Create unit test when it is necessary. A problem is only considered resolved when the related unit test passes.

 
