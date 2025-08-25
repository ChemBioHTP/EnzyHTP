# Plan

The task is to make the `_parameterize_modified_res` more robust and use accurate force field.
1. Before searh the ncaa library, we want to check if this maa is already supported by the force field user specified. Make a seperate function for this and call it. 
    - This new function will be `check_residue_name_ff_support(res_code: str, force_fields: List[str], work_dir: Path) -> bool`
    - res_code is the 3-letter code of the residue.
    - In this function, we use a simple tleap call to verify the support. 
        ```
        # source all user specified force fields
        source leaprc.protein.ff14SB
        ...

        # put the res code in 
        model = sequence { RES }

        # print out the description
        desc model

        quit
        ```
        explore what will happen with supported/unsupported case.
    - make unit test for this function and do test driven development
    - put it in `_parameterize_modified_res` and configure all the arguments using information we know in `_parameterize_modified_res`
2. It is possible that user will use force field combinations such as ff19SB and ff19SB_modAA. Make sure `get_protein_force_field` if given such list find the correct canonincal AA protein force field. Add this case in the unit test and address any problem.

