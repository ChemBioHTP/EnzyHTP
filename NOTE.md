I wish to implement umbrella sampling in EnzyHTP. Can you help me implement this feature?

You should write **concise**, **modular** code that extends and **follows all the conventions of EnzyHTP**. Please look carefully at the codebase before starting to implement umbrella sampling. Importantly, your addition should be consistent with existing functions in EnzyHTP.

In particular, consider the science API `equi_md_sampling` and how it is currently used. The umbrella sampling science API should be similar to existing science APIs.

I have written for you some scripts that already can perform umbrella sampling. However, they need to be modularized.

experiments/exp-umbrella/run18/template_md_only.py
experiments/exp-umbrella/run18/umbrella_figures.py

## Features to add

The following should be created:
- a class that generalizes reaction coordinate + constraint + value called CollectiveVariable at `...`
- a file that interfaces WHAM at `...`
- a file that handles umbrella sampling at `...`
- a file that handles umbrella sampling analysis at `...`

The following requirements should be met at minimum:
1. A function to perform umbrella sampling, given a reaction coordinate. Parameters should be tunable, like spring constant.
2. A function to calculate the Probability Density with respect to reaction coordinate (umbrella_probability_density.png)
3. A function to calculate the Potential of Mean Force via WHAM. 
4. WHAM should be integrated as an external dependency with its own interface. (this is an external dependency, I know that EnzyHTP handles those systematically but I'm not sure exactly how, can you investigate?)
5. A function to plot the PMF with respect to reaction coordinate.
6. A function to inspect phase space to ensure that there is phase space overlap.

## Completion Condition

1. All features defined above are implemented.
2. Unit tests are written that covers every feature you add. 
3. Integration tests are written that covers every feature you add.
4. Tests that consider edge cases such as:
    - what if phase space overlap is insufficient?
    - ...
5. A test that reproduces the results from the scripts that I provided you. 
6. All tests should pass.
