### Final Revised Plan

**Task:** Refactor prepin to mol2 conversion to preserve atom types. We need to replace the current way for converting a prepin file to a mol2 file (line 406) because it will mess up the atom type.

**Approach:** A new `MolDescData` data class will be introduced to serve as a standardized, intermediate representation for molecular description data, ensuring the I/O modules remain decoupled.

1.  **Define a new `MolDescData` class:**
    *   A new file will be created: `enzy_htp/structure/mol_desc_data.py`.
    *   This file will contain a `dataclass` named `MolDescData`.
    *   This class will store generic molecular information, including:
        *   `name: str`
        *   `atoms: List[Dict]` (with keys: `id`, `atom_name`, `atom_type`, `charge`, `coords`)
        *   `bonds: List[Dict]` (with keys: `atom1_id`, `atom2_id`, `bond_type`)

2.  **Create a `prepin` to `MolDescData` converter:**
    *   In `enzy_htp/structure/structure_io/prepin_io.py`, a new method `get_mol_desc_data(path: str) -> MolDescData` will be added to `PrepinParser`.
    *   Use `PrepinParser._parse_prepin_file` and translate the data into a `MolDescData` object.

3.  **Create a `MolDescData` to `mol2` writer:**
    *   In `enzy_htp/structure/structure_io/mol2_io.py`, a new method `write_from_mol_desc_data(data: MolDescData, outfile: str)` will be added to `Mol2Parser`.
    *   This method will take a `MolDescData` object, generate the corresponding `.mol2` file content, and write it to `outfile`.

4.  **Implement the orchestrator method in `AmberInterface`:**
    *   In `enzy_htp/_interface/amber_interface.py`, a new method `convert_mol_desc_format(self, old_file: str, new_file: str, old_format: str, new_format: str)` will be created.
    *   This method will act as a dispatcher for file format conversions.
    *   Initially, it will only support `prepin` to `mol2`. It will check the `old_format` and `new_format` arguments and raise a `NotImplementedError` for any other combination.
    *   It will first call `PrepinParser.get_mol_desc_data(old_file)` to get the intermediate `MolDescData` object.
    *   Then, it will call `Mol2Parser.write_from_mol_desc_data(data, new_file)` to write the final `.mol2` file.

5.  **Update `AmberParameterizer` to use the new method:**
    *   In `enzy_htp/_interface/amber_interface.py`, inside the `_parameterize_modified_res` method of the `AmberParameterizer` class, the call to `self.parent_interface.run_antechamber` for the prepin-to-mol2 conversion (around line 406) will be replaced.
    *   The new call will be `self.parent_interface.convert_mol_desc_format(old_file=prepin_path, new_file=mol_desc_path, old_format='prepin', new_format='mol2')`.
