# Feedback to Claude Code

1. The dsi unit test should have at least one unit test that does not mock the underlying engine.
2. The `_residue_list_to_amber_mask`, should have an argument to let user choose if hydrogens are excluded.
3. The `_residue_list_to_amber_mask` should call and reuse `get_interval_str_from_list` instead of write long and hard to read loops that re-implement similar function
4. Amber seems never allow using different chain ids. So when the same residue number with different chain ids shows up in `_residue_list_to_amber_mask` give a warning.
5. In `_residue_list_to_amber_mask`, when there are multiple residues, it should be ":10,15,20&!@H=" instead of ":10,:15,:20&!@H=". 
6. In `calculate_dsi_metrics`, the ensemble will not always contain Amber format `topology_source_file` and `coordinate_list`, the current solution is to use
    ```
    self.convert_top_to_prmtop(stru_esm.topology_source_file, tmp_prmtop_path)
    self.convert_traj_to_nc(stru_esm.coordinate_list, tmp_nc_path, topology_path=tmp_prmtop_path)
    ```
    See how they are used in `get_coord_covariance`
7. We dont need `test_calculate_dsi_metrics_mock`, make it `test_calculate_dsi_metrics` and don't mock anything there. See how `test_get_coord_covariance` create the ensemble for testing.
