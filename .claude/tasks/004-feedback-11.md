# Feedback for Claude Code

1. In `af2_predict`, we want to provide more comprehensive output from the af2 run, and make it the return value of the function. Here is the detailed design:
   {
    <sequence/tuple_of_sequence(from the input)> : {
        "best_model" : Structure,
        "best_model_plddt" : [float, ...],
        "best_model_pae" : [float, ...],
        "best_model_index" : int,
        "model_1" : Structure,
        "model_1_plddt" : [float, ...],
        "model_1_pae" : [float, ...],
        "model_2" : Structure,
        "model_2_plddt" : [float, ...],
        "model_2_pae" : [float, ...],
        ...
    },
    ...
   }
Use `test_af2_predict_colabfold` to verify your change.
