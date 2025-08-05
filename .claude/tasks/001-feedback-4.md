# Feedback for Claude Code

Several problems in your code.
    1. The run prepgen could hide error in the try block. I have edited the code to fix it. Review it and make sure it works. (i.e.: add a unit test in the case that will fail)
    2. The existing unit test of run_prepgen need a rework. It now leaves intermediate files in the working dir and will not correct raise a failed prepgen run.
