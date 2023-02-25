'''
Misc helper func and class
'''
from distutils.command.config import config
import math
from subprocess import CompletedProcess, SubprocessError, run
import time
import os
from typing import List
import numpy as np

from Class_Conf import Config
'''
====
Tree
====
'''
class Child():
    def __init__(self):
        self.parent = None
    def set_parent(self, parent_obj):
        self.parent = parent_obj

        return self

'''
Text
'''
line_feed = Config.line_feed

'''
file system
'''
def mkdir(dir):
    if os.path.exists(dir):
        pass
    else:
        os.makedirs(dir)

def is_empty_dir(dir_path):
    '''
    check if the dir_path is an empty dir
    '''
    if os.path.isdir(dir_path):
        if not os.listdir(dir_path):
            return 1
        else:
            return 0
    else:
        print(f"No such directory: {dir_path}")
        return 2

'''
math
'''
def set_distance(p1,p2,d):
    '''
    return a coord of p3
    -- p3 --
    origin:     p1
    direction: (p1,p2)
    distance:   d
    '''
    p1 = np.array(p1)
    p2 = np.array(p2)
    #direction vector
    v1 = (p2 - p1)/np.linalg.norm(p1 - p2)
    p3 = p1 + v1 * d

    return tuple(p3)


def get_distance(p1,p2):
    d1 = np.array(p2) - np.array(p1)
    D = np.linalg.norm(d1)

    return D


def get_field_strength_value(p0, c0, p1, p2=None, d1=None):
    '''
    return field strength E of *p0(c0)* at *p1* in direction of *p2-p1* or *d1*
    -- E = kq/r^2 -- (Unit: kcal/mol)
    point charge:   c0 in p0 
    point:          p1
    direction:      p2-p1 or d1
    '''
    # Unit
    k = 332.4   # kcal*Ang/(mol*e^2) = (10^10)*(4.184^-1)*(Na)*(10^-3)*(1.602*10^-19)^2 * 9.0*10^-9 N*m^2/C^2
    q = c0                      # e
    p0 = np.array(p0)           # Ang 
    p1 = np.array(p1)           # Ang
    if d1 == None:
        d1 = np.array(p2) - p1
    else:
        d1 = np.array(d1)
    d1 = d1/np.linalg.norm(d1)  # Ang

    # Get r
    r = p1 - p0
    r_m = np.linalg.norm(r)
    r_u = r/r_m

    # Get E
    E = (k * q / (r_m**2)) * r_u
    # Get E in direction
    Ed = np.dot(E, d1)          # kcal/(mol*e*Ang)

    return Ed


def get_center(p1, p2):
    '''
    return the center of p1 and p2
    '''
    p3 = 0.5 * (np.array(p1) + np.array(p2))

    return tuple(p3)

def round_by(num: float, cutnum: float) -> int:
    '''
    round the float number up if the decimal part is larger than cutnum
    otherwise round down
    '''
    dec_part, int_part = math.modf(num)
    if dec_part > cutnum:
        int_part += 1
    return int(int_part)

'''
Cheminfo
'''
def Conformer_Gen_wRDKit(input_mol, out_path, numConfs=50):
    '''
    Generate small molecular conformers using RDKit. (wrote in 2022/1/13 for ReactiveDocking)
    If input 3D structure, require protonated one. (treat add H in a different function)
    ---
    In: SMILES | MOL2 | PDB | SDF
    Out: SDF | PDB
    ---
    numConfs: number of conformers output
    '''
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from Class_Conf import Config

    # input type support 
    # [input_mol --> mol]
    # file or SMILES
    if not os.path.exists(input_mol):
        if Config.debug > 1:
            print('WARNING: Conformer_Gen_wRDKit: input is not a file. Using SMILES mode.')
        # SMILES
        mol = Chem.MolFromSmiles(input_mol)
        mol = Chem.AddHs(mol)
    else:
        # file
        sfx = input_mol.split('.')[-1].strip()
        if sfx == 'smiles':
            mol = Chem.MolFromSmiles(input_mol)
            mol = Chem.AddHs(mol)
        if sfx == 'mol2':
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
    o_sfx = out_path.split('.')[-1].strip()
    if o_sfx == 'sdf':
        with Chem.SDWriter(out_path) as of:
            Chem.rdmolfiles.SDWriter.SetKekulize(of, False)
            for i in cids:
                of.write(mol, confId=i)
    if o_sfx == 'pdb':
        with Chem.PDBWriter(out_path) as of:
            for i in cids:
                of.write(mol, confId=i)
        
    return out_path


'''
Rosetta
'''
def generate_Rosetta_params(input_file, out_dir, resn, out_pdb_name, if_conformer=0, overwrite=0):
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
'''
def decode_atom_mask(stru, mask, ifsolvent=0):
    '''
    decode atom mask and return a list of atom ids.
    Base on the correponding structure obj.
    ===Only support residues in whole now & do not discrimminate chains start id count from 1===
    '''
    atom_ids = []

    resi_ids = mask[1:].strip().split(',')
    for i in range(len(resi_ids)-1, -1, -1):
        if '-' in resi_ids[i]:
            r1 = int(resi_ids[i].split('-')[0])
            r2 = int(resi_ids[i].split('-')[1])
            resi_ids.extend(list(range(r1, r2+1)))
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
    '''
    use repr() to store data
    expect a eval() to decode the stored file
    '''
    tag = repr(tag)
    
    with open(out_path, 'a') as of:
        of.write('===TAG==='+line_feed)
        of.write(tag+line_feed)
        for item in data.keys():
            of.write('---'+item+'---'+line_feed)
            of.write(repr(data[item])+line_feed)

    return out_path

def chunked(iter_, size):
    '''
    chunk iter_ by size and return generator for chunked list 
    '''
    return (iter_[position : position + size] for position in range(0, len(iter_), size))

def get_localtime(time_stamp=None):
    if time_stamp is None:
        return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    else:
        return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time_stamp))

def run_cmd(cmd, try_time=1, wait_time=3, timeout=120) -> CompletedProcess:
    '''
    try running the info cmd {try_time} times and wait {wait_time} between each run if subprocessexceptions are raised.
    default be 1 run.
    along with common run() settings (including exception handling)
    # TODO(shaoqz): should use this as general function to run commands in local shell.
    '''
    for i in range(try_time):
        try:
            this_run = run(cmd, timeout=timeout, check=True,  text=True, shell=True, capture_output=True)
        except SubprocessError as e:
            if Config.debug > 0:
                print(f'Error running {cmd}: {repr(e)}')
                print(f'    stderr: {str(e.stderr).strip()}')
                print(f'    stdout: {str(e.stdout).strip()}')
                print(f'trying again... ({i+1}/{try_time})')
        else: # untill there's no error
            if Config.debug > 0:
                if i > 0:
                    print(f'finished {cmd} after {i+1} tries @{get_localtime()}')
            return this_run
        # wait before next try
        time.sleep(wait_time)
    # exceed the try time
    raise SubprocessError(f'Failed running `{cmd}` after {try_time} tries @{get_localtime()}')
    # TODO change to a custom error

def delete_idx_line(target_line_str: str, idx: int) -> str:
    """split the string by line feed and deleted the indexed line and put
    them back to a string"""
    lines = target_line_str.split(os.linesep)
    lines.pop(idx)
    return "\n".join(lines)

def check_complete_metric_run(mutant: List[str], data_file_path: str) -> bool:
    """Check if a mutation is fully finished in a typical enzy_htp run"""
    data_dict_list = extract_enzy_htp_data(data_file_path)
    for data_dict in data_dict_list:
        muta_flags = ["".join(x) for x in data_dict["TAG"]]
        if mutant == muta_flags:
            return True
    return False

def extract_enzy_htp_data(data_file_path: str) -> List[dict]:
    """extract typical enzy_htp run data that generated by write_data()"""
    result = []
    with open(data_file_path) as f:
        lines = f.readlines()
        d_flag = 0
        m_flag = 0
        for i, line in enumerate(lines):
            if i == 0:
                m_flag = 1
                m_data = {}
                continue
            if i+1 == len(lines):
                if d_flag:
                    line_data = eval(line.strip())
                    m_data[Term] = line_data
                    d_flag = 0
                result.append(m_data)
                break
            if 'TAG' in line:
                result.append(m_data)
                m_flag = 1
                m_data = {}
                continue
            if m_flag:
                MutaFlag = eval(line.strip())
                m_data['TAG'] = MutaFlag
                m_flag = 0
                continue

            if line.strip()[:3]+line.strip()[-3:] == '------':
                Term=line.strip().strip('-')
                d_flag = 1
                continue
            if d_flag:
                line_data = eval(line.strip())
                m_data[Term] = line_data
                d_flag = 0
                continue
    return result
