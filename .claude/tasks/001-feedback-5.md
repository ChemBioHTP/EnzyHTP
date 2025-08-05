# Feedback for Claude Code

Several problems in your code.
    1. `test_run_prepgen` need a rework. You need to generate a .mc file and keep it in the data dir. (i.e. using a seperate run.) And avoid generating the .mc file during the test. I have deleted everything in that unit test. Please plan what to test first and complete it and make it work.
    2. For fixing the problem isolated in `test_structure_deepcopy_isolation`. Please 
       a. read and understand the current custom deepcopy method Structure class use
       b. read and understand how connectivity is stored in Structure
       c. based on your understanding analyze the reason of the problem
       d. write your analysis of the reason and a concrete plan about what changes should be made on the custom deepcopy function under `.claude/tasks/001-claude-plan.md`
