"""TODO(CJ)

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-06-11
"""
# TODO(CJ): make the unit test for this part
import re
from typing import List, Dict, Set
from enzy_htp.structure import Structure, Residue


def decode_mask_amber(mask: str) -> List[int]:
    # TODO(CJ): documentation
    mask = mask.strip()
    if not mask.startswith(":"):
        raise TypeError()

    result = mask[1:].strip().split(",")
    result = list(map(float, result))
    return result


def create_selection(struct: Structure, atom_mask: str, fix_end: str = None):
    # def get_sele_list(self, atom_mask, cap_strat='H', prepi_path=None):
    """
    interface with class ONIOM_Frame. Generate a list for sele build. Make sure use same pdb as the one generate the frame.
    ------------
    resi_list: selected residue list
    atom_mask: atom selection with the standard grammer of Amber (incomplete)
    cap_strat: fix valence of the cut bond. (default: H)
            - H: add H to where the original connecting atom is.
                special fix for classical case:
                "sele by residue" (cut N-C) -- adjust dihedral for added H on N.
            = Interface with write_sele_lines:
                add {fix_flag+element_mark: coord} in sele_lines
    ------------
    return a sele list:
    - Fixing atoms are labeled as qm_atom_id-qm_atom_cnt_id-distance
    - backbone atoms are marked as b at the end (for qmcluster charge calculation)
    - other atoms use _ as place holder
    return a sele map:
    (PDB atom id -> QM atom id)
    """
    fix_end = "H"
    sele_lines = {}
    # decode atom_mask (maybe in helper later) TODO
    resi_list = atom_mask[1:].strip().split(",")
    all_resi_list = struct.residues

    # decode and get obj
    sele_stru_objs = []
    for resi in resi_list:
        chain_id = re.match("[A-Z]", resi)
        resi_id = int(re.match("[0-9]+", resi).group(0))
        if chain_id == None:
            for resi in all_resi_list:
                if resi_id == resi.idx():
                    resi_obj = resi
        else:
            chain_id = chain_id.group(0)
            resi_obj = self.chains[int(chain_id) - 65]._find_resi_id(resi_id)

        sele_stru_objs.append(resi_obj)

    # combine the sele
    sele_atoms = []
    for obj in sele_stru_objs:
        for atom in obj.atoms():
            sele_atoms.append(atom)

    print(sele_atoms)
    if fix_end != None:
        self.get_connect(prepi_path=prepi_path)
    # operate on the sele objs
    for atom in sele_atoms:
        # add current atom
        atom.get_ele()
        if type(atom.parent) != Ligand:
            if atom.name in ["C", "CA", "O", "N", "H", "HA"]:
                sele_lines[str(atom.id) + "b"] = atom.ele
            else:
                sele_lines[str(atom.id) + "_"] = atom.ele
        else:
            sele_lines[str(atom.id) + "_"] = atom.ele

        if fix_end != None:
            # search for cut bond
            for cnt_atom in atom.connect:
                if not cnt_atom in sele_atoms:
                    if fix_end == "H":
                        d_XH = X_H_bond_length[atom.name]
                        label = "-".join(
                            (str(atom.id), str(cnt_atom.id), str(d_XH)))
                        fix_atom = "H"
                    if fix_end == "Me":
                        # TODO
                        pass
                    # write to sele_lines
                    sele_lines[label] = fix_atom
    # make sele_map (PDB atom id -> QM atom id)
    sele_map = {}
    for i, key in enumerate(sele_lines.keys()):
        if key[-1] not in "1234567890":
            key = key[:-1]
        sele_map[key] = i + 1

    if Config.debug >= 1:
        print("Selected QM cluster atoms: ")
        print(sele_lines)

    return sele_lines, sele_map
