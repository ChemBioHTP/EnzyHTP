"""Aids in the preparation of enzyme systems by providing an initial, unrefined "seed" for a given ligand. Use the below
functions:
    
    + seed_with_coordinates(): 
    + seed_with_transplants():
    + seed_with_constraints():
    + seed_with_analog():
    + seed_with_pdb_structure();

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2024-05-27
"""
import json
from typing import List, Tuple, Dict, Set, Callable, Union

from rdkit import Chem as _rchem
import numpy as np
import pandas as pd
import enzy_htp as eh
from pathlib import Path
from copy import deepcopy
from collections import defaultdict

from enzy_htp import interface, config
from enzy_htp.core import file_system as fs

from enzy_htp.core import _LOGGER

import enzy_htp.chemical as chem
from enzy_htp.structure import PDBParser, Structure, Ligand, Mol2Parser, Chain, Atom
from enzy_htp.structure.structure_constraint import StructureConstraint
from enzy_htp._interface import Mole2Cavity

from enzy_htp.structure.structure_operation import (
    atom_name_similarity
)

from .ligand_moves import (
    mimic_torsions,
    mimic_torsions_mcs,
    ligand_mcs_score
)

def seed_with_coordinates(ligand:Ligand,
                            coords:Union[Tuple[float, float, float], List[Tuple[float,float,float]]],
                            rng_seed:int=1996,
                            work_dir:str=None) -> None:
    """Seeds the location of the Ligand()'s center of mass using the supplied seeds. If multiple seeds are supplied, one
    is chosen at random. 

    Args:
        ligand: The Ligand() to seed.
        coords: A (float, float, float) or List[(float, float, float)] to choose the seed from.
        rng_seed: Random number generation seed. Optional.
        work_dir: temp directory where the work is done.  Optional.


    Returns:
        Nothing.
    """
    np.random.seed( rng_seed )

    if not isinstance( coords, list ):
        coords = [ coords ]

    coord = np.array(np.random.choice( coords ))
    current = np.array(ligand.geom_center)
    ligand.shift(
        coord - current
    )


def seed_with_transplants(ligand:Ligand,
                     similarity_metric:str,
                     similarity_cutoff:float = 0.45,
                     use_cache:bool=False,
                     work_dir:str=None
) -> Ligand:
    """Implementation function that places a specified ligand into a structure using the alphafill 
    algorithm (https://doi.org/10.1038/s41592-022-01685-y). 

    Args:
        ligand: The Ligand() object to seed.
        simliarity_metric: The method by which similarity is scored. Valid options include: atom_names or mcs
        similaritY_cutoff: The cutoff value for the similarity metric. Default is 0.45.
        use_cache: Should a file cache be used to avoid calling alphafill again? Default is False.
        work_dir: Location where temporary and output files should be saved. Optional

    Returns:
        An aligned, deepcopied Ligand().
    """
    _LOGGER.info(f"Beginnning placement of ligand {ligand} into the structure...")
    _LOGGER.info(f"Calling out to AlphaFill to fill structure...")
    # atom_names mcs
    
    (filled_structure, df) = interface.alphafill.fill_structure(
        ligand.parent.parent,use_cache=use_cache, work_dir=work_dir)
    _LOGGER.info("Filled structure using AlphaFill!")
    
    ligand_atoms:Set[str]=set([aa.name for aa in ligand.atoms])
    
    similarity_score:List[float] = list()
    sim_func = None
    if similarity_metric == 'atom_names':
        sim_func = atom_name_similarity
    elif similarity_metric == 'mcs':
        sim_func = ligand_mcs_score
    else:
        err_str:str=f"The similarity metric {similarity_metric} is not supported!"
        _LOGGER.error(err_str)
        raise TypeError(err_str)
    
    df['similarity_score'] = df.apply(lambda row: sim_func( ligand, row.mol ), axis=1 )
    similarity_mask = df.similarity_score >= similarity_cutoff

    if not similarity_mask.sum():
        _LOGGER.error(f"No transplants exist with similarity_cutoff of {similarity_cutoff:.3f}! Exiting...")
        exit( 1 )
    
    df = df[similarity_mask].reset_index(drop=True)
    similarity = df.similarity_score.to_numpy()
    df = df[np.isclose(similarity, np.max(similarity))].reset_index(drop=True)
    df.sort_values(by='clash_count',inplace=True)
    df.reset_index(drop=True, inplace=True)
    template:Ligand = df.iloc[0].mol
    
    if not ligand.is_ligand():
        ligand.atoms[0].coord = template.atoms[0].coord
    else:
        if df.iloc[0].similarity_score >= 0.95:
            interface.rdkit.apply_coords(template, ligand) 
        elif similarity_metric == 'atom_names':
            mimic_torsions( template, ligand )
        elif similarity_metric == 'mcs':
            mimic_torsions_mcs( template, ligand )

def seed_with_analog(ligand:Ligand,
                analog_template:Ligand,
                analog:Ligand,
                seed_atom:Atom,
                analog_atom:Atom,
                work_dir:str=None) -> None:

    """Seed a Ligand using an Analog Ligand. Aligns the seed Atom from that Ligand to the specified analog Atom.

    Args:
        ligand: The Ligand to be seeded.
        analog_template: The analog Ligand that will serve as a template
        analog: The analog Ligand whose coordinates will be used.
        seed_atom: The seed Atom that will be aligned to the analog atom.
        analog_atom: The analog Atom whose coordinates will be used to place the seed Atom.
        work_dir: Directory where all work will be done. Optional.

    Returns:
        Nothing.
    """

    if not work_dir:
        work_dir = config['system.SCRATCH_DIR']


    mimic_torsions(analog_template, analog)

    target_location = None

    for atom in analog.atoms:
        if atom.name == analog_atom.name:
            target_location = np.array(atom.coord)
            break
    
    assert target_location is not None

    current_location = None
    for atom in ligand.atoms:
        if atom.name == seed_atom.name:
            current_location = np.array(atom.coord)
            break

    assert current_location is not None
    shift = target_location - current_location

    ligand.shift( shift )

def seed_with_constraints(ligand:Ligand,
        constraints:List[StructureConstraint], 
        delta:float=0.75,
        work_dir:str=None ) -> None:
    """Seed a Ligand using StructureConstraints. First uses Mole2 to identify all cavity spaces, 
    then loops through a 3D grid of that space, choosing the seed location that doesn't clash with existing Atom()'s
    and has the best constraint score, i.e. that best satisfies the constraints.

    Args:
        ligand: The Ligand() to seed. 
        constraints: The StructureConstraint()'s to be used for scoring purposes.
        delta: The delta to create the point grid, used for x/y/z axes.
        work_dir: Where the temporary files should be saved. Optional.

    Returns:
        Nothing.

    """
    stru = ligand.parent.parent
    assert constraints
    if not work_dir:
        work_dir = config['system.SCRATCH_DIR']

    relevant_constraints: List[str] = list()

    cavities:List[Mole2Cavity] = interface.mole2.identify_cavities(stru, work_dir=work_dir)

    seed_locations:List = list()
    for cc in cavities:
        vert_matrix = np.array(cc.points())
        (x_min, y_min, z_min) = np.min(vert_matrix,axis=0)
        (x_max, y_max, z_max) = np.max(vert_matrix,axis=0)
        contained_points = list()
        x_vals = list(np.arange(x_min, x_max+delta, delta))
        y_vals = list(np.arange(y_min, y_max+delta, delta))
        z_vals = list(np.arange(z_min, z_max+delta, delta))
        
        candidates = list()
        for x in x_vals:
            for y in y_vals:
                for z in z_vals:
                    candidates.append(np.array([x, y, z]))
        candidates = np.array(candidates)
        for included, cp in zip(cc.contains_points(candidates), candidates):
            if included:
                seed_locations.append(cp)


    seed_locations = np.array(seed_locations)
    scores = list()
    
    for sl in seed_locations:
    
        lig_start = ligand.geom_center
        shift = sl - lig_start

        ligand.shift(shift)

        curr_energy = 0.0
        for cst in constraints:
            if sum([aa.parent == ligand for aa in cst.atoms]) > 1:
                continue
            curr_energy += interface.rosetta.score_energy(cst)
        scores.append( curr_energy )            

    scores = np.array(scores)

    score_mask = np.isclose(scores, np.min(scores))

    clash_counts = list()
    
    for slidx,sl in enumerate(seed_locations[score_mask]):
        lig_start = ligand.geom_center
        shift = sl - lig_start

        ligand.shift(shift)
        clash_ct = 0

        for res in stru.residues:
            if res == ligand:
                continue
            clash_ct += res.clash_count( ligand )
    
        clash_counts.append( clash_ct ) 

    clash_counts = np.array(clash_counts)

    clash_mask = np.isclose(clash_counts, np.min(clash_counts))

    final_locations = seed_locations[score_mask][clash_mask]
    seed = final_locations[0]

    ligand.shift( seed - ligand.geom_center)


def seed_with_pdb_structure(ligand:Ligand,
                     pdb_code:str,
                     ligand_sele:str,
                     align_sele:str=None,
                     work_dir:str=None
) -> None:
    """Seed a Ligand using coordinates from an existing protein databank database structure. 
    Requires at least one heavy atom name match between candidate and template Ligand()'s.

    Args:
        ligand: The Ligand to be seeded.
        pdb_code: The four-letter PDB code of the reference entry. Case-insensitive.
        ligand_sele: The pymol-formatted selection for the reference ligand in the refence PDB entry.
        align_sele: The pymol-formatted selection that the Ligand()'s parent structure will be aligned to.
        work_dir: Directory where the work will be done. Optional.

    Returns:
        Nothing.
    """
    _LOGGER.info(f"Beginnning placement of ligand {ligand} into the structure...")
    _LOGGER.info(f"Calling out to AlphaFill to fill structure...")
    # atom_names mcs
    session = interface.pymol.new_session() 
    stru = ligand.parent.parent
    (stru_sele, session) = interface.pymol.load_enzy_htp_stru(session, stru)
   
    if align_sele is None:
        align_sele = ''
    else:
        align_sele = f"and ({align_sele})"
        
    temp_file = "./temp.mol2"
    args = [
        ('fetch', pdb_code),
        ('remove', 'solvent'),
        ('align', f"{pdb_code} {align_sele}", f"{stru_sele} {align_sele}"),
        ('save', temp_file, ligand_sele)
        ]
    interface.pymol.general_cmd(session, args )
    lp = Mol2Parser()
    template = lp.get_ligand( temp_file )
    
    ligand_atoms:Set[str]=set([aa.name for aa in ligand.atoms])
    template_atoms:Set[str]=set([aa.name for aa in template.atoms]) 
    ligand_atoms_no_h:Set[str]=set([aa.name for aa in ligand.atoms if aa.element !='H'])
    template_atoms_no_h:Set[str]=set([aa.name for aa in template.atoms if aa.element !='H']) 
    
    if ligand_atoms_no_h.issubset(template_atoms_no_h) or template_atoms_no_h.issubset(ligand_atoms_no_h):
        from rdkit.Chem import AllChem
        template_to_ligand = dict()
        atom_names = list()
        tmol = interface.rdkit.mol_from_ligand(template, False)
        lmol = interface.rdkit.mol_from_ligand(ligand, False)
        for aidx, at in enumerate(tmol.GetAtoms()):
            target_name:str=at.GetPropsAsDict()['_TriposAtomName']
            atom_names.append(target_name)
            for lidx, al in enumerate(lmol.GetAtoms()):
                if target_name == al.GetPropsAsDict()['_TriposAtomName']:
                    template_to_ligand[aidx] = lidx
                    break
            else:
                pass
        
        lconf = lmol.GetConformer()
        tconf = tmol.GetConformer()
        
        for tidx, lidx in template_to_ligand.items():
            
            lconf.SetAtomPosition(lidx, tconf.GetAtomPosition(tidx))

        interface.rdkit.update_ligand_positions(ligand, lmol)
    else:
        err_str="seed_with_pdb_structure(): No heavy atom match between candidate and template Ligand()'s" 
        _LOGGER.error(err_str)
        raise TypeError(err_str)

