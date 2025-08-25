# Feedback to Claude Code

1) In `get_residue_pka_from_stru` 
   1. the work_dir should be by default the SCRATCH dir of enzy_htp. You can reference how other module use it. (same for `get_residue_pka_from_pdb`)
   2. Let's put the logic about remove solvent in residue_pka.py so that there are no circular import nor in-function import.

2) In `get_residue_pka_from_pdb`
   1. Change the name of this function to run_propka so that it is a main wrapper of propka.
   2. In this case, we should extract the full output information for each residue, including these my other project uses:
        ```
        row_dict["res_num"] = atom.res_num
        row_dict["ins_code"] = atom.icode
        row_dict["res_name"] = atom.res_name
        row_dict["chain_id"] = atom.chain_id
        row_dict["group_label"] = group.label
        row_dict["group_type"] = getattr(group, "type", None)
        row_dict["pKa"] = group.pka_value
        row_dict["model_pKa"] = group.model_pka
        row_dict["buried"] = group.buried
        if group.coupled_titrating_group:
            row_dict["coupled_group"] = group.coupled_titrating_group.label
        else:
            row_dict["coupled_group"] = None
        ```
    
3) In `_filter_target_residues`
   1. The current way of just using res_num to align target residues to the result from propka is not robust. It is possible that two residues have the same residue number but on different chains. In this case, we need chain id to align them as well.
   2. As 1 need chain id, we should no longer omit them in the saved structure, we will just rely on removing solvent to avoid the chain id overflow problem.

