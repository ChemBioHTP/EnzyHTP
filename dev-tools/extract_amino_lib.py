"""a tool for extracting information from Amber amino_xxx.lib file

Usage:
1. write the main()
2. use argv"""
import sys
from typing import Dict, List, Any

def read_from_lib_file(target_file: str) -> List[Dict]:
    """read atomic charge from the target lib file"""
    with open(target_file) as f:
        entry_list = f.read().split('!')
        entry_list = entry_list[3:]
        result = []

        for entry in entry_list:
            entry_lines = entry.strip().split('\n')
            entry_data_type = entry_lines[0].split(' ')[1]
            if entry_data_type == "table":
                entry_labels = entry_lines[0].split('  ') # entry.ALA.unit.positions table  dbl x  dbl y  dbl z
                entry_ids = entry_labels[0].split(' ')[0]
                entry_labels = entry_labels[1:]
            else:
                entry_labels = None
                entry_ids = entry_lines[0].split(' ')[0]
            entry_ids = entry_ids.split(".") # entry.ALA.unit.positions
            entry_resi = entry_ids[1]
            entry_type = entry_ids[-1]
            if len(entry_resi) == 4: # C/N-ter
                entry_resi = entry_resi[1:]
            entry_data = {
                "residue": entry_resi,
                "type": entry_type,
                "values": {
                    "type" : entry_data_type,
                    "labels" : entry_labels,
                    "value_lines" : entry_lines[1:],
                }
            }
            result.append(entry_data)
    return result

def decode_type_atoms(entry: Dict) -> List[Dict]:
    """decode for entry type: atoms
    name  type  typex  resx  flags  seq  elmnt  chg"""
    result = []
    entry = entry["values"]
    labels = [lb.split(" ")[1] for lb in entry["labels"]]
    for line in entry["value_lines"]:
        line_parts = line.strip().split(" ")
        line_dict = {}
        for k, v in zip(labels, line_parts):
            line_dict[k] = v
        result.append(line_dict)
    return result
        
def fetch_resi_atomic_charge(target_file: str) -> Dict[str, Dict[str, float]]:
    """workflow for fetch atomic charge information form the target file"""
    result = {}
    entry_list = read_from_lib_file(target_file)
    entry_list = filter(lambda x: x["type"] == "atoms", entry_list)
    for entry in entry_list:
        resname = entry["residue"]
        atoms_data = decode_type_atoms(entry)
        res_charge_map = {}
        for atom in atoms_data:
            res_charge_map[atom["name"].strip('"')] = float(atom["chg"])
        result[resname] = res_charge_map
    return result

def write_result(obj: Any, out_file: str):
    """write result object to temp file for copying"""
    with open(out_file, "w") as of:
        of.write(repr(obj))

def main():
    write_result(fetch_resi_atomic_charge(sys.argv[1]), sys.argv[2])

if __name__ == "__main__":
    main()
