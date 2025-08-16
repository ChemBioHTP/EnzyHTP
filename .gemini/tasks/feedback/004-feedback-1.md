# Feedback for Gemini

0. I made some changes to the code you wrote. Be aware of it.
1. The prediction.py need a dictionary to contain all supported engines. See analysis/dsi.py for example.
2. For run() in alphafold_interface actually implement it. We need to reference analysis/binding.py and support the user choice of 1) running alphafold locally or 2) submit it as a HPC job as a ClusterJob.
3. For implementing the run() method, reference the following scratch code, fit them in EnzyHTP's design principle:
    ```

    ```
4. You also need unit tests for alphafold_interface.
