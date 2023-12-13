import numpy as np
from scipy.spatial.transform import Rotation as R
from typing import List, Tuple, Dict
from ..structure import Structure, Residue, Atom

#TODO(CJ): add methylamide -> cterm, acetate -> nterm

def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

def get_h_atom(end:str) -> List[Atom]:
    if end == 'nterm':
        aname = 'HP11'
    elif end == 'cterm':
        aname = 'HP21'

    return [
            Atom({'atom_name': aname, 'x_coord':0.000, 'y_coord': 0.000, 'z_coord':  0.000, 'atom_number':1})
    ]

class _ResidueDummy:
    def __init__(self, name:str):
        self.name = name

def get_ch3_atoms(end:str, res_dummy:str) -> List[Atom]:
    """TODO(CJ)"""
    if end == 'nterm':
        anames = "CP1 HP11 HP12 HP13".split()
    elif end == 'cterm':
        anames = "CP2 HP21 HP22 HP23".split()

    return [
            Atom({'atom_name': anames[0], 'x_coord':0.000, 'y_coord': 0.000, 'z_coord':  0.000,'atom_number':1},parent=res_dummy),
            Atom({'atom_name': anames[1], 'x_coord':0.360, 'y_coord':-1.029, 'z_coord':  0.000,'atom_number':1},parent=res_dummy),
            Atom({'atom_name': anames[2], 'x_coord':0.360, 'y_coord': 0.514, 'z_coord':  0.891,'atom_number':1},parent=res_dummy),
            Atom({'atom_name': anames[3], 'x_coord':0.360, 'y_coord': 0.514, 'z_coord': -0.891,'atom_number':1},parent=res_dummy)
    ]

def needs_nterm_capping(res:Residue, stru:Structure, residue_list:List[Tuple[str,int]]) -> bool:
    tup_key:Tuple[str,int] = (res.parent.name, res.idx)
    return stru.has_residue( f"{tup_key[0]}.{tup_key[1]-1}") and (tup_key[0], tup_key[1]-1) not in residue_list

def needs_cterm_capping(res:Residue, stru:Structure, residue_list:List[Tuple[str,int]]) -> bool:
    tup_key:Tuple[str,int] = (res.parent.name, res.idx)
    return stru.has_residue( f"{tup_key[0]}.{tup_key[1]+1}") and (tup_key[0], tup_key[1]+1) not in residue_list

def cap_residue(res:Residue, stru:Structure, strategy:str, side:str) -> List[Atom]:
    """TODO(CJ)"""
   
    chain1, num1, = res.key()
    if side == 'nterm':
        atom_key1 = f"{chain1}.{num1}.N"
        atom_key2 = f"{chain1}.{num1-1}.C"
    elif side == 'cterm':
        atom_key1 = f"{chain1}.{num1}.C"
        atom_key2 = f"{chain1}.{num1+1}.N"

    start_atom:Atom = stru.get_atom(atom_key1)
    end_atom:Atom = stru.get_atom(atom_key2)

    p1 = np.array(start_atom.coord)
    p2 = np.array(end_atom.coord)
    
    direction = p2 - p1
    direction /= np.linalg.norm(direction)
    
    cap:List[Atom] = None
    if strategy == 'CH3':
        cap = get_ch3_atoms(side, res)
    elif strategy == 'H':
        cap = get_h_atom(side)
   
    bond_distance:float = get_bond_distance(strategy, side)

    rot_mat = rotation_matrix_from_vectors(np.array([1,0,0]), direction)
    
    for aa in cap:
        pt = np.array(aa.coord)
        pt = np.transpose(np.matmul(rot_mat, np.transpose( pt  )))
        pt += (direction*bond_distance) + p1 
        aa.coord = pt

    return cap

def get_bond_distance(strategy:str, side:str) -> float:
    """ """
    #TODO(CJ): update this: need methylamide and acetyl
    _BOND_DISTANCE_MAPPER:Dict[Tuple[str,str], float] = {
        ('H', 'nterm'):1.1,
        ('H', 'cterm'):1.0,
        ('CH3', 'nterm'):1.5,
        ('CH3', 'cterm'):1.5,
    }


    result:float = _BOND_DISTANCE_MAPPER.get((strategy, side), None)

    if result is None:
        _LOGGER.error(f"No bond distance for {strategy}-capping strategy on {side} side. Supported methods include:")
        for (strat, side_name) in _BOND_DISTANCE_MAPPER.keys():
            _LOGGER.error(f"strategy: {strat}, side: {side_name}")
        _LOGGER.error("Exiting...")
        exit( 1 )

    return result

