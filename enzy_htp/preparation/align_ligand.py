"""Functionality for aligning the coordinates of a ligand to a target ligand using backbone dihedrals. Only the fuction
align_ligand() should be called. All other functions are for implementation purposes only and SHOULD NOT be called
directly by the user.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-10-09
"""

import pandas as pd
import numpy as np
import rdkit.Chem as _rchem
from rdkit.Chem import rdMolTransforms as rdmt
from collections import defaultdict
from typing import List, Dict, Set, Tuple
from pathlib import Path

from enzy_htp import config, interface

from enzy_htp.core import file_system as fs
from enzy_htp.core import _LOGGER 


def align_ligand(template:str, ligand:str, work_dir:str=None, outfile:str=None) -> str:
    """Method that aligns a ligand to a template geometry using backbone dihedrals. Relies on the atom naming scheme
    used in the PDB. This approach enables greater flexibility as ATP, ADP, and AMP can all be aligned to a template
    ADP molecule. Derivatized molecules can also be aligned to a target substrate.

    Args:
        template: str() filename containing the target molecule geometry.
        ligand: str() filename containing the original molecule.
        work_dir: the directory to do work in. Optional. Defaults to config['system.SCRATCH_DIR'].
        outfile: The path to save the ligand to. Defaults to <template>_aligned.mol2

    Returns:
        The path to the aligned ligand.
    """
    if work_dir is None:
        work_dir = config['system.WORK_DIR']

    t_info:Dict = _molecule_info(template, template=True)
    l_info:Dict = _molecule_info(ligand, template=False)

    _adjust_torsions(t_info, l_info, adjust_rings=True)
    _adjust_torsions(t_info, l_info, adjust_rings=False)

    session = interface.pymol.new_session()
    args = [('load', ligand)]
    for aidx, atom in enumerate(l_info['mol'].GetAtoms()):
        pos = l_info['mol'].GetConformer().GetAtomPosition(aidx)

        args.append(
            ('alter_state', -1, f'name {l_info["mapper"][aidx][aidx]}', f'(x,y,z) = ({pos.x},{pos.y},{pos.z})')
        )
    
    temp_path=Path(ligand)
    if outfile is None:
        outfile:str=f"{work_dir}/{temp_path.stem}_placed{temp_path.suffix}"
    
    obj1:str=f"{Path(template).stem}"
    obj2:str=temp_path.stem

    atoms:str=" or ".join(map(lambda nn: f"name {nn}",  set(t_info['names']).intersection(set(l_info['names']))))

    args.extend([
        ('load', template),
        ('align', f"{obj2} and ({atoms})", f"{obj1} and ({atoms})"),        
        ('save', outfile, obj2)
    ])

    interface.pymol.general_cmd(session, args)

    return outfile


def non_h_bonds(atom : _rchem.Atom) -> int:
    """How many bonds does the molecule have where the other atom is NOT Hydrogen?"""
    count = 0
    for bb in atom.GetBonds():
        if bb.GetBeginAtom().GetAtomicNum() == 1 or bb.GetEndAtom().GetAtomicNum() == 1:
            pass
        else:
            count += 1
    return count            

def _map_atoms(mol:_rchem.RWMol, fname;str, cutoff:float=0.001) -> Tuple[Dict,List[str]]:
    """Function that maps atom indices from rdkit to PDB-style atom names. Returns a dict()
    containing the individual mappings as well as a list() of the PDB-style atom names as str()'s.

    Args:
        mol: The molecule to map to.
        fname: The str() name of the file the molecule came from. 
        cutoff: The float() distance allowed for direct atom mapping. Defaults to 0.001.

    Returns:
        A Tuple() with the layout (mapper, atom names).

    """
    session = interface.pymol.new_session()
    interface.pymol.general_cmd(session, [('delete', 'all'),('load', fname)])
    df:pd.DataFrame=interface.pymol.collect(session, 'memory', 'name x y z'.split())
    mapper = dict()        
    points = np.transpose(np.array([df.x.to_numpy(), df.y.to_numpy(), df.z.to_numpy()]))
    names:List[str] = df.name.to_list()
    
    terminals = set()      
    
    for aidx, atom in enumerate(mol.GetAtoms()):
        temp = {'terminal':False}
        pos = mol.GetConformer().GetAtomPosition(aidx)
        point = np.array([pos.x, pos.y, pos.z])
        if non_h_bonds(atom) == 1:
            temp['terminal'] = True
        dists = np.sqrt(np.sum((points-point)**2,axis=1))
        it = np.argmin(dists)
        if dists[it] >= cutoff:
            continue
        temp[atom.GetIdx()] = names[it]
        temp[names[it]] = atom.GetIdx()

        temp['bond_ct'] = len(atom.GetBonds())
        
        mapper[atom.GetIdx()] = temp
        mapper[names[it]] = temp

    return (mapper,names)

def _get_non_h_bonds(parent,partner):
    bonds = list()
    for bb in parent.GetBonds():
        a1, a2 = bb.GetBeginAtom(), bb.GetEndAtom()

        if a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1:
            continue

        if a1.GetIdx() == parent.GetIdx() and a2.GetIdx() == partner.GetIdx():
            continue
        
        if a2.GetIdx() == parent.GetIdx() and a1.GetIdx() == partner.GetIdx():
            continue
   
        if bb.GetIsAromatic():
            continue

        bonds.append(bb)
    
    return bonds 

def _select_bond(anchor, bonds, mapper) -> _rchem.Bond:
    """

    Args:
        anchor:
        bonds:
        mapper:

    Returns:
        The selected rdkit.Chem.Bond() object.
    """

    if len(bonds) == 1:
        return bonds[0]

    counts = list()
    atype = list()
    for bb in bonds:
        if bb.GetBeginAtom().GetAtomicNum() == 1 or bb.GetBeginAtom().GetAtomicNum() == 1:
            continue
        other = bb.GetOtherAtom(anchor)
        o_idx = other.GetIdx()
        counts.append(mapper[o_idx]['bond_ct'])
        atype.append(other.GetAtomicNum())

    counts = np.array(counts)
    atype = np.array(atype)

    max_ct = np.max(counts)

    if sum(counts==max_ct)==1:
        return bonds[np.argmax(counts)]
    
    temp = np.array(bonds)[counts==max_ct]
    atype = atype[counts==max_ct]
    
    max_atype=np.max(atype)

    if sum(max_atype==atype)==1:
        return temp[np.argmax(max_atype)]

    return temp[0]


def _ring_neighbors(mol,atoms):
    bond_mapper=defaultdict(list)

    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetAtomicNum() == 1 or  bond.GetEndAtom().GetAtomicNum() == 1:
            continue

        bond_mapper[bond.GetBeginAtom().GetIdx()].append(bond.GetEndAtom().GetIdx())
        bond_mapper[bond.GetEndAtom().GetIdx()].append(bond.GetBeginAtom().GetIdx())
    
    rings = list()
    for pp in _rchem.GetSymmSSSR(mol):
        rings.append(set(list(pp)))

    cut_points=list()
    for aidx in [atoms[1],atoms[2],atoms[0],atoms[3]]:
        ring = None
        for rr in rings:
            if aidx in rr:
                ring = rr                

        if not ring:
            continue

        for idx in bond_mapper[aidx]:
            if idx in atoms:
                continue
            if idx not in ring:
                continue
            
            cut_points.append((aidx,idx))
                

    return cut_points

def _get_backbone_torsions(mol:_rchem.RWMol, mapper:Dict) -> pd.DataFrame:
    """Gets a torsions DataFrame for ONLY the backbone bonds in the supplied molecule. Similar
    to _get_all_torsions() but only gets dihedral atom names for backbone bonds.

    Args:
        mol: The rdkit.Chem.RWMol() to get dihedral names for.
        mapper: An atom mapper from the _map_atoms() function.

    Returns:
        A DataFrame containing backbone dihedral atom names as well as ring proximity information.
    """
    dihedrals:List[List[int]] = list()
    mapped_dihedrals:List[List[str]] = list()
    bond_idxs = list()
    in_ring = list()
    ring_neighbors = list()
    
    for bond in mol.GetBonds():

        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()

        if (a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1):
            continue

        a1_mapper, a2_mapper = mapper[a1.GetIdx()], mapper[a2.GetIdx()]
        
        if (a1_mapper['terminal'] or a2_mapper['terminal']):
            continue
        
        a1_bonds=_get_non_h_bonds(a1,a2)
        a2_bonds=_get_non_h_bonds(a2,a1)

        if (not a1_bonds) or (not a2_bonds):
            continue

        selected1, selected2 = _select_bond(a1,a1_bonds,mapper), _select_bond(a2,a2_bonds,mapper)

        a0=selected1.GetOtherAtom(a1)        
        a3=selected2.GetOtherAtom(a2)
        atoms=[a0.GetIdx(),a1.GetIdx(),a2.GetIdx(),a3.GetIdx()]
        mapped_names=list(map(lambda aa: mapper[aa][aa], atoms))

        if mapped_names[0] < mapped_names[-1]:
            mapped_names = mapped_names[::-1]
            atoms = atoms[::-1]

        in_ring.append(a1.IsInRing() and a2.IsInRing() and not (
            a0.GetIsAromatic() and a1.GetIsAromatic() and a2.GetIsAromatic() and a3.GetIsAromatic()
        ))

        neighbors = list()
        if in_ring[-1]:
            neighbors = _ring_neighbors(mol, atoms)

        ring_neighbors.append(neighbors)

        dihedrals.append(tuple(atoms))
        mapped_dihedrals.append(tuple(mapped_names))


    return pd.DataFrame({
            'dihedrals':dihedrals,
            'mapped_dihedrals':mapped_dihedrals,
            'in_ring':in_ring,
            'ring_neighbors':ring_neighbors
            })

def _get_all_torsions(mol:_rchem.RWMol, mapper:Dict) -> pd.DataFrame:
    """Gets a torsions DataFrame for ALL the bonds in the supplied molecule. Similar
    to _get_backbone_torsions() but gets ALL possible dihedral angle names..

    Args:
        mol: The rdkit.Chem.RWMol() to get dihedral names for.
        mapper: An atom mapper from the _map_atoms() function.

    Returns:
        A DataFrame containing all dihedral atom names as well as ring proximity information.
    """

    dihedrals:List[List[int]] = list()
    mapped_dihedrals:List[List[str]] = list()
    bond_idxs = list()
    in_ring = list()
    ring_neighbors = list()
    
    for bond in mol.GetBonds():

        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()

        if (a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1):
            continue

        a1_mapper, a2_mapper = mapper[a1.GetIdx()], mapper[a2.GetIdx()]
        
        if (a1_mapper['terminal'] or a2_mapper['terminal']):
            continue
        
        a1_bonds=_get_non_h_bonds(a1, a2)
        a2_bonds=_get_non_h_bonds(a2, a1)

        if (not a1_bonds) or (not a2_bonds):
            continue

        selected1, selected2 = _select_bond(a1, a1_bonds, mapper), _select_bond(a2, a2_bonds, mapper)

        for b1 in a1_bonds:
            for b2 in a2_bonds:
                
                a0=b1.GetOtherAtom(a1)        
                a3=b2.GetOtherAtom(a2)
                atoms:List[int]=[a0.GetIdx(),a1.GetIdx(),a2.GetIdx(),a3.GetIdx()]
                mapped_names=list(map(lambda aa: mapper[aa][aa], atoms))
                
                if mapped_names[0] < mapped_names[-1]:
                    mapped_names = mapped_names[::-1]
                    atoms = atoms[::-1]

                in_ring.append(a1.IsInRing() and a2.IsInRing() and not (
                    a0.GetIsAromatic() and a1.GetIsAromatic() and a2.GetIsAromatic() and a3.GetIsAromatic()
                ))
                
                neighbors = list()
                if in_ring[-1]:
                    neighbors = _ring_neighbors(mol, atoms)
        
                ring_neighbors.append(neighbors)
        
                dihedrals.append(tuple(atoms))
                mapped_dihedrals.append(tuple(mapped_names))

    return pd.DataFrame({
            'dihedrals':dihedrals,
            'mapped_dihedrals':mapped_dihedrals,
            'in_ring':in_ring,
            'ring_neighbors':ring_neighbors
            })


def _molecule_info(fname:str, template:bool=True) -> Dict:  
    """Gathers information about the molecule needed for alignment. Information includes the dihedral angles,
    atom names, and other information. If the molecule is specified as a template, only backbone dihedral angles
    are taken. Otherwise all possible dihedrals are found. The returned dict() has the following keys with 
    the respective information:
        
        + 'df': DataFrame containing each dihedral as well as information on ring proximity.
        + 'mol': The rdkit.Chem.RWMol() molecule object for the file.
        + 'mapper': A dict() which maps atom indices or names to a dict() which contains atom bonding information.
        + 'names': A list() of atom names for the molecule.

    Args:
        fname: Filename of the molecule to analyze as a str().
        template: Is this a template ligand?
    
    Returns:
        A dict() with the keys 'df', 'mol', 'mapper', and 'names'.

    """
    mol:_rchem.Molecule = interface.rdkit._load_molecule(fname, sanitize=True)
    (mapper, names) = _map_atoms(mol, fname)
    if template:
        df:pd.DataFrame = _get_backbone_torsions(mol,mapper) 
    else:
        df:pd.DataFrame = _get_all_torsions(mol,mapper) 

    return {
        'df':df,
        'mol':_rchem.RWMol(mol),
        'mapper':mapper,
        'names':names
    }


def _adjust_torsions(template:Dict, ligand:Dict, adjust_rings=False) -> None:
    """
    Args:

    Returns:
        Nothing.
    """
    #print(template)
    intersection:Set[str] = set(template['df'].mapped_dihedrals.to_list()).intersection(set(ligand['df'].mapped_dihedrals.to_list()))

    tconf=template['mol'].GetConformers()[0]
    lconf=ligand['mol'].GetConformers()[0]

    for dihedral in sorted(list(intersection)):
        t_row = template['df'][template['df'].mapped_dihedrals==dihedral].iloc[0]
        l_row = ligand['df'][ligand['df'].mapped_dihedrals==dihedral].iloc[0]
        
        if l_row.in_ring != adjust_rings:
            continue

        n1,n2=None,None
        if l_row.in_ring:   
            for pr in l_row.ring_neighbors:
                if pr is not None:
                    (n1,n2)=pr
                    break
            
            if n1 and n2:
                ligand['mol'].RemoveBond(n1,n2)
                _rchem.SanitizeMol(ligand['mol'])
                lconf=ligand['mol'].GetConformers()[0]
        
        rdmt.SetDihedralDeg(lconf,*l_row.dihedrals,
            rdmt.GetDihedralDeg(tconf,*t_row.dihedrals)
        )

        if n1 and n2:
            ligand['mol'].AddBond(n2,n1)
            _rchem.SanitizeMol(ligand['mol'])
            lconf = ligand['mol'].GetConformers()[0]
            _rchem.AllChem.MMFFOptimizeMolecule(ligand['mol'])

