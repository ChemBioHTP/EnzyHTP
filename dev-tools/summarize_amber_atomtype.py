
import pandas as pd
from parmed import load_file
import os
import pprint
from collections import Counter

# --- Step 1: Build the detailed master map ---

mapping_path = "/panfs/accrepfs.vampire/home/shaoq1/temp/amberff_lib_dat_mapping.csv"
df_mapping = pd.read_csv(mapping_path)
main_chain_atoms = ["N", "H", "CA", "C", "O", "OXT"]
master_mapping = {}

for index, row in df_mapping.iterrows():
    leaprc_file = row["leaprc File"]
    lib_files_str = row["lib Files"]
    if not isinstance(lib_files_str, str): continue
    lib_files = [f.strip() for f in lib_files_str.split(",")]
    master_mapping[leaprc_file] = {}
    for lib_file in lib_files:
        full_path = f"/gpfs51/dors2/csb/apps/amber22/dat/leap/lib/{lib_file}"
        if not os.path.exists(full_path): continue
        try:
            lib = load_file(full_path)
            for res_name, residue in lib.items():
                if res_name not in master_mapping[leaprc_file]:
                    master_mapping[leaprc_file][res_name] = {}
                for atom in residue.atoms:
                    if atom.name in main_chain_atoms:
                        master_mapping[leaprc_file][res_name][atom.name] = atom.type
        except Exception:
            continue
# --- Step 2: Analyze the master map to create the summarized dictionary ---

FORCEFIELD_SUMMARY = {}

# Define standard amino acids (non-terminal) to determine the mainstream mapping
standard_amino_acids = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'CYX', 'GLU', 'GLN', 'GLY', 'HIS', 
    'HID', 'HIE', 'HIP', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
    'THR', 'TRP', 'TYR', 'VAL'
}

def get_residue_description(res_names):
    res_names = set(res_names)
    if not res_names: return "Unknown"

    # Check for specific capping groups first
    if res_names == {'ACE'}: return "N-terminal Acetyl cap"
    if res_names == {'NME'}: return "C-terminal N-Methylamide cap"

    # Check for terminal residue patterns, ensuring the base residue is a standard one
    is_c_term = all(r.startswith('C') and len(r) == 4 and r[1:] in standard_amino_acids for r in res_names)
    if is_c_term: return "C-terminal residues"

    is_n_term = all(r.startswith('N') and len(r) == 4 and r[1:] in standard_amino_acids for r in res_names)
    if is_n_term: return "N-terminal residues"
    
    # Check for proline-related, which might have unique backbone patterns
    if all('PRO' in r for r in res_names): return "Proline-related residues"
    
    # Default to listing the residues
    residue_list_str = ", ".join(sorted(list(res_names)))
    return f"Specific residues: {residue_list_str}"

for ff, res_maps in master_mapping.items():
    if not res_maps: continue

    # --- Define Mainstream Mapping ---
    # Count occurrences of each full backbone mapping within standard amino acids.
    standard_res_maps_counter = Counter(
        frozenset(res_map.items()) for res_name, res_map in res_maps.items()
        if res_name in standard_amino_acids and len(res_map) >= 5 # Ensure it's a full backbone map
    )
    
    mainstream_mapping = {}
    if standard_res_maps_counter:
        # The mainstream mapping is the most common one found.
        most_common_mapping_fs = standard_res_maps_counter.most_common(1)[0][0]
        mainstream_mapping = dict(most_common_mapping_fs)

    # --- Classify all residues based on the mainstream mapping ---
    mainstream_residues = []
    non_mainstream_groups = {}

    for res_name, res_map in res_maps.items():
        if not res_map: continue # Skip empty entries
            
        # A residue is mainstream if its mapping is identical to the mainstream one.
        if res_map == mainstream_mapping:
            mainstream_residues.append(res_name)
        # Otherwise, it's non-mainstream. Group it with others that have the same mapping.
        else:
            mapping_key = frozenset(res_map.items())
            if mapping_key not in non_mainstream_groups:
                non_mainstream_groups[mapping_key] = []
            non_mainstream_groups[mapping_key].append(res_name)

    # --- Format the non-mainstream summary ---
    non_mainstream_summary = []
    for mapping_fs, res_list in non_mainstream_groups.items():
        description = get_residue_description(res_list)
        non_mainstream_summary.append({
            'description': description,
            'residues': sorted(res_list),
            'mapping': dict(mapping_fs)
        })

    # --- Assemble the final summary for the force field ---
    FORCEFIELD_SUMMARY[ff] = {
        'mainstream_mapping': mainstream_mapping,
        'mainstream_residues': sorted(mainstream_residues),
        'non_mainstream_mappings': sorted(non_mainstream_summary, key=lambda x: x['description'])
    }

# --- Step 3: Write the final dictionary to a .py file ---

output_string = f"FORCEFIELD_SUMMARY = {pprint.pformat(FORCEFIELD_SUMMARY)}\n"
output_path = "/panfs/accrepfs.vampire/home/shaoq1/temp/amberff_summary.py"
with open(output_path, "w") as f:
    f.write(output_string)

print(f"Final summary dictionary has been written to {output_path}")
