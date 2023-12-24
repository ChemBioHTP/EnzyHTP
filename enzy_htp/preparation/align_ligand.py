"""Functionality for aligning the coordinates of a ligand to a target structure using mapped PDB-style atom names. Note that atoms 
without names will not be easily aligned. Only align_ligand() should be called. All other functions are for implementation purposes 
only and SHOULD NOT be called directly by the user.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-10-09
"""

from collections import defaultdict
from typing import List, Dict, Set, Tuple
from pathlib import Path

import pandas as pd
import numpy as np
import rdkit.Chem as _rchem
from rdkit.Chem import rdMolTransforms as rdmt

from enzy_htp import config, interface

from enzy_htp.core import file_system as fs
from enzy_htp.core import _LOGGER


def align_ligand(template: str,
                 ligand: str,
                 work_dir: str = None,
                 outfile: str = None,
                 minimize_its: int = 10000,
                 mapping_tol: float = 0.001,
                 cst_tol: float = 0.05,
                 cst_penalty: float = 100.0,
                 bond_tol: float = 0.05) -> str:
    """Method that aligns a ligand to a template geometry using atomic locations. Relies on the atom naming scheme
    used in the PDB. This approach enables greater flexibility as ATP, ADP, and AMP can all be aligned to a template
    ADP molecule, for example. Derivatized molecules can also be aligned to a target substrate. The general workflow
    consists of the below steps:
        1. mapping of atoms between template structure and ligand structure
        2. moving ligand structure coordinates to template structure locations
        3. setting distorted bond distances back to original distance
        4. constrainted minimization of non-mapped atoms

    Args:
        template: str() filename containing the target molecule geometry. Must be in .mol2 format.
        ligand: str() filename containing the original molecule. Must be in .mol2 format.
        work_dir: the directory to do work in. Optional. Defaults to config['system.SCRATCH_DIR'].
        outfile: The path to save the ligand to. Defaults to <work_dir>/<template>_aligned.mol2
        minimize_its: Number of iterations allowed in constrained optimization. Defaults to 10000.
        mapping_tol: Tolerance in Angstroms allowed when mapping atoms between template and ligand. Defaults to 0.001
        cst_tol: Allowed positional distance during minimzation for constrained atoms. Defaults to 0.05.
        cst_penalty: Energetic penalty for positional distance deviations during minimization for constrained atoms. Defaults to 100.0.
        bond_tol: Tolerance in bond distances caused by mapping that triggers manipulation. Defaults to 0.05. 

    Returns:
        The path to the aligned ligand.
    """
    err_msg: List[str] = list()
    if fs.get_file_ext(template) != '.mol2':
        err_msg.append(f"\tThe supplied file {template} is not a .mol2 file!")

    if fs.get_file_ext(ligand) != '.mol2':
        err_msg.append(f"\tThe supplied file {template} is not a .mol2 file!")

    if err_msg:
        _LOGGER.error(f"Errors in input files:")
        for em in err_msg:
            _LOGGER.error(em)
        _LOGGER.error("Exiting...")
        exit(1)

    t_info: Dict = _molecule_info(template, mapping_tol)
    l_info: Dict = _molecule_info(ligand, mapping_tol)

    _rchem.SanitizeMol(l_info['mol'])
    props = _rchem.AllChem.MMFFGetMoleculeProperties(l_info['mol'])
    ff = _rchem.AllChem.MMFFGetMoleculeForceField(l_info['mol'], props)

    orig_mat = _rchem.rdmolops.Get3DDistanceMatrix(l_info['mol'])

    mapped_atoms = set()
    for bb in t_info['bonded'].intersection(l_info['bonded']):
        mapped_atoms.add(bb[0])
        mapped_atoms.add(bb[1])

    t_conf = t_info['mol'].GetConformers()[0]
    l_conf = l_info['mol'].GetConformers()[0]

    for an in mapped_atoms:
        l_conf.SetAtomPosition(l_info['mapper'][an], t_conf.GetAtomPosition(t_info['mapper'][an]))
        ff.MMFFAddPositionConstraint(l_info['mapper'][an], cst_tol, cst_penalty)

    updated_mat = _rchem.rdmolops.Get3DDistanceMatrix(l_info['mol'])

    for bond in l_info['bonded']:
        bond_obj=l_info['bond_mapper'][bond]
        if bond_obj.IsInRing():
            continue
        idx1, idx2 = l_info['mapper'][bond[0]], l_info['mapper'][bond[1]]
        avg = (updated_mat[idx1][idx2] + orig_mat[idx1][idx2]) / 2
        if abs(updated_mat[idx1][idx2] - orig_mat[idx1][idx2]) / avg >= bond_tol:
            if bond[0] in mapped_atoms:
                i, j = idx1, idx2
            else:
                j, i = idx1, idx2
            rdmt.SetBondLength(l_conf, i, j, orig_mat[idx1][idx2])

    ff.Initialize()
    ff.Minimize(maxIts=minimize_its)
    temp_path = Path(ligand)

    if work_dir is None:
        work_dir = config['system.WORK_DIR']

    if outfile is None:
        outfile: str = f"{work_dir}/{temp_path.stem}_placed{temp_path.suffix}"

    session = interface.pymol.new_session()
    args = [('load', ligand)]
    for aidx, atom in enumerate(l_info['mol'].GetAtoms()):
        pos = l_info['mol'].GetConformer().GetAtomPosition(aidx)

        args.append(('alter_state', -1, f'name {l_info["mapper"][aidx]}', f'(x,y,z) = ({pos.x},{pos.y},{pos.z})'))

    args.append(('save', outfile))

    interface.pymol.general_cmd(session, args)

    return outfile


def _map_atoms(mol: _rchem.RWMol, fname: str, cutoff: float = 0.001) -> Tuple[Dict, List[str]]:
    """Function that maps atom indices from rdkit to PDB-style atom names. Returns a dict()
    containing the individual mappings. Atoms are mapped in both directions, i.e. name->index and index->name.

    Args:
        mol: The molecule to map to.
        fname: The str() name of the file the molecule came from. 
        cutoff: The float() distance allowed for direct atom mapping. Defaults to 0.001.

    Returns:
        The mapper dict() which maps string atom names to atom index and vice-versa.
    """
    session = interface.pymol.new_session()
    interface.pymol.general_cmd(session, [('delete', 'all'), ('load', fname)])
    df: pd.DataFrame = interface.pymol.collect(session, 'memory', 'name x y z'.split())
    points = np.transpose(np.array([df.x.to_numpy(), df.y.to_numpy(), df.z.to_numpy()]))
    names: List[str] = df.name.to_list()

    mapper = dict()

    for aidx, atom in enumerate(mol.GetAtoms()):
        pos = mol.GetConformer().GetAtomPosition(aidx)
        point = np.array([pos.x, pos.y, pos.z])
        dists = np.sqrt(np.sum((points - point)**2, axis=1))
        it = np.argmin(dists)

        if dists[it] >= cutoff:
            continue
        mapper[atom.GetIdx()] = names[it]
        mapper[names[it]] = atom.GetIdx()

    return mapper


def _molecule_info(fname: str, dist_tol: float) -> Dict:
    """Gathers information about the molecule needed for alignment. Information includes the dihedral angles,
    atom names, and other information. If the molecule is specified as a template, only backbone dihedral angles
    are taken. Otherwise all possible dihedrals are found. The returned dict() has the following keys with 
    the respective information:
        
        + 'mol': The rdkit.Chem.RWMol() molecule object for the file.
        + 'mapper': A dict() which maps atom indices or names to a dict() which contains atom bonding information.
        + 'bonded': A set() containing all mapped atom names as tuples with format (aname_1, aname_2).

    Args:
        fname: Filename of the molecule to analyze as a str().
        dist_tol: The mapping distance tolerance in Angstroms.

    
    Returns:
        A dict() with the keys 'mol', 'mapper', and 'bonded'. Descriptions of values are given above.
    """
    mol: _rchem.Molecule = interface.rdkit._load_molecule(fname, sanitize=True)
    if mol is None:
        mol = interface.rdkit._load_molecule(fname)

    mapper = _map_atoms(mol, fname, dist_tol)

    bonded = set()

    bond_mapper = dict()

    for bond in mol.GetBonds():
        a1_idx, a2_idx = bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()
        bonded.add((mapper[a1_idx], mapper[a2_idx]))
        bonded.add((mapper[a2_idx], mapper[a1_idx]))
        bond_mapper[(mapper[a1_idx], mapper[a2_idx])] = bond
        bond_mapper[(mapper[a2_idx], mapper[a1_idx])] = bond

    return {
        'mol': _rchem.RWMol(mol),
        'mapper': mapper,
        'bonded': bonded,
        'bond_mapper': bond_mapper
    }
