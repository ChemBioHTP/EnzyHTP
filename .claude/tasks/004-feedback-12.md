# Feedback for Claude Code

1. I have prepared an example output dir from the native AF2 install @test/_interface/data/test_af2_output_ref/ . Create a unit test based on this example dir and see if the current `_parse_comprehensive_results` works well with it. If not, revise the function to support it while keeping the colabfold support. (We no longer want pae in the output; plddt is in the .pkl file under the 'plddt' key)
