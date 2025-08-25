# Feedback for Claude Code

Several problems in your code.
    1. The `run_prepgen` could hide error in the try block. I have edited the code to fix it. Review it and make sure it works. (i.e.: add a unit test in the case that will fail)
    2. The existing unit test of run_prepgen need a rework. It now leaves intermediate files in the working dir and will not correct raise a failed prepgen run.
    3. `test_structure_deepcopy_isolation` should not be in `test_amber_interface` because it has nothing to do with Amber. Put it in the right place.
    4. Propose a plan to fix the deepcopy problem isolated in `test_structure_deepcopy_isolation`. Write the plan in a new file under `.claude/tasks/001-claude-plan-{index}.md`. You should come up with several plans index them in the filename.
