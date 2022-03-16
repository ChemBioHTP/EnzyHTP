"""
Misc helper func and class
"""
from distutils.command.config import config
from AmberMaps import Resi_list
import os
import numpy as np
"""
====
Tree
====
"""


class Child:

    def __init__(self):
        self.parent = None

    def set_parent(self, parent_obj):
        self.parent = parent_obj

        return self


"""
Text
"""
line_feed = "\n"
"""
func
"""


def mkdir(dir):
    if os.path.exists(dir):
        pass
    else:
        os.makedirs(dir)


"""
math
"""


def set_distance(p1, p2, d):
    """
    return a coord of p3
    -- p3 --
    origin:     p1
    direction: (p1,p2)
    distance:   d
    """
    p1 = np.array(p1)
    p2 = np.array(p2)
    # direction vector
    v1 = (p2 - p1) / np.linalg.norm(p1 - p2)
    p3 = p1 + v1 * d

    return tuple(p3)


def get_distance(p1, p2):
    d1 = np.array(p2) - np.array(p1)
    D = np.linalg.norm(d1)

    return D

    #---
    #numConfs: number of conformers output
    '''
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from Class_Conf import Config

    # input type support
    # [input_mol --> mol]
    # file or SMILES
    if not os.path.exists(input_mol):
        if Config.debug > 1:
            print(
                "WARNING: Conformer_Gen_wRDKit: input is not a file. Using SMILES mode."
            )
        # SMILES
        mol = Chem.MolFromSmiles(input_mol)
        mol = Chem.AddHs(mol)
    else:
        # file
        sfx = input_mol.split(".")[-1].strip()
        if sfx == "smiles":
            mol = Chem.MolFromSmiles(input_mol)
            mol = Chem.AddHs(mol)
        if sfx == "mol2":
            mol = Chem.MolFromMol2File(input_mol, removeHs=0)
            Chem.rdmolops.AssignAtomChiralTagsFromStructure(mol)
        if sfx == 'pdb':
            mol = Chem.MolFromPDBFile(input_mol, removeHs=0)
            Chem.rdmolops.AssignAtomChiralTagsFromStructure(mol)
        if sfx == 'sdf':
            mol = Chem.SDMolSupplier(input_mol, removeHs=0)[0]
            Chem.rdmolops.AssignAtomChiralTagsFromStructure(mol)
        
    # calculate conformers & minimize
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=numConfs, numThreads=0, pruneRmsThresh=0.1)
    AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0)

    # output
    o_sfx = out_path.split(".")[-1].strip()
    if o_sfx == "sdf":
        with Chem.SDWriter(out_path) as of:
            Chem.rdmolfiles.SDWriter.SetKekulize(of, False)
            for i in cids:
                of.write(mol, confId=i)
    if o_sfx == "pdb":
        with Chem.PDBWriter(out_path) as of:
            for i in cids:
                of.write(mol, confId=i)

    return out_path


#Rosetta
'''


def generate_Rosetta_params(input_file,
                            out_dir,
                            resn,
                            out_pdb_name,
                            if_conformer=0,
                            overwrite=0):
    '''
    generate rosetta params file using molfile_to_params.py
    The command will be:
    [some_path]/molfile_to_params.py -n [resn] -p [out_pdb_name] --keep-names
    This command will output files under the working directory. The function will move them to out_dir
    ------
    resn:           if resn is 'same' then extract the resi name from the input_file
    out_pdb_name:   if out_pdb_name is 'same' then extract the resi name from the input_file
    if_conformer:   if add --conformers-in-one-file flag
    overwrite:      if add --clobber flag
    '''
    pass
    # return params_path, out_pdb_path, conformers_path (can only under the same dir as params)


'''
misc
"""


def decode_atom_mask(stru, mask, ifsolvent=0):
    """
    decode atom mask and return a list of atom ids.
    Base on the correponding structure obj.
    ===Only support residues in whole now & do not discrimminate chains start id count from 1===
    """
    atom_ids = []

    resi_ids = mask[1:].strip().split(",")
    for i in range(len(resi_ids) - 1, -1, -1):
        if "-" in resi_ids[i]:
            r1 = int(resi_ids[i].split("-")[0])
            r2 = int(resi_ids[i].split("-")[1])
            resi_ids.extend(list(range(r1, r2 + 1)))
            del resi_ids[i]
    # clean
    resi_ids = [int(i) for i in resi_ids]
    resi_ids.sort()

    all_resi = stru.get_all_residue_unit(ifsolvent=ifsolvent)
    for r_id in resi_ids:
        for resi in all_resi:
            if resi.id == r_id:
                for atom in resi:
                    atom_ids.append(atom.id)

    return atom_ids


def write_data(tag, data, out_path):
    """
    use repr() to store data
    expect a eval() to decode the stored file
    """
    tag = repr(tag)

    with open(out_path, "a") as of:
        of.write("===TAG===" + line_feed)
        of.write(tag + line_feed)
        for item in data.keys():
            of.write("---" + item + "---" + line_feed)
            of.write(repr(data[item]) + line_feed)

    return out_path
        
Mutant assigner
'''
#TODO
