"""Aids in the preparation of enzyme systems by providing an initial, unrefined "seed" for a given ligand. Use the below
functions:
    
    + seed_with_coordinates(): 
    + seed_with_transplants():

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
from enzy_htp.structure import PDBParser, Structure, Ligand, Mol2Parser, Chain
from enzy_htp.structure.structure_constraint import StructureConstraint
from enzy_htp._interface import Mole2Cavity

from enzy_htp.structure.structure_operation import (
    ligand_mcs,
    atom_name_similarity
)

from .ligand_moves import mimic_torsions
from enzy_htp.geometry import minimize as geo_minimize


def seed_with_coordinates(ligand:Ligand,
                            coords:Union[Tuple[float, float, float], List[Tuple[float,float,float]]],
                            rng_seed:int=1996,
                            minimize:bool=False,
                            min_iter:int=1,
                            work_dir:str=None) -> None:
    """Seeds the location of the Ligand()'s center of mass using the supplied seeds. If multiple seeds are supplied, one
    is chosen at random. Ligand() can be minimized at end.

    Args:
        ligand: The Ligand() to seed.
        coords: A (float, float, float) or List[(float, float, float)] to choose the seed from.
        rng_seed: Random number generation seed. Optional.
        minimize: Should Ligand()-only minimization be performed?
        min_iter: How many minimization iterations should be performed? Optional.
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

    if minimize:
        minimize_ligand_only( ligand, min_iter, work_dir )


def seed_ligand(stru: Structure, 
    ligand:Ligand, 
    method:str, 
    minimize:bool=False, 
    constraints:List[StructureConstraint]=None, 
    work_dir:str=None, 
    **kwargs) -> None:
    """

    Args:
        ligand:
        method:
        minimize:
        constraints:
        work_dir:

    Returns:
        Nothing.
    """

    if not work_dir:
        work_dir = config['system.SCRATCH_DIR']

    fs.safe_mkdir(work_dir)
    
    func:Callable = SEEDING_METHOD_MAPPER.get( method )

    if func is None:
        err_msg:str=f"The supplied method '{method}' is not supported. Supported include: {', '.join(SEEDING_METHOD_MAPPER.keys())}"
        _LOGGER.error(err_msg)
        raise ValueError(err_msg)

    kwargs['constraints'] = constraints
    func(stru, ligand, work_dir=work_dir, **kwargs)

    if minimize:
        sele = f'(chain {ligand.parent.name} and resi {ligand.idx})'
        #sele = f'byres polymer.protein within 2 of {lig_sele} or {lig_sele}'
        geo_minimize(stru, movemap=[
            {'sele':sele, 'bb':'true', 'chi':'true', 'bondangle':'true'}, 
            {'sele':f'not ({sele})', 'bb':'false', 'chi':'false', 'bondangle':'false'}, 

        ],n_iter=1)



def seed_with_transplants(ligand:Ligand,
                     similarity_metric:str,
                     similarity_cutoff:float = 0.45,
                     minimize:bool=False,
                     min_iter:int=1,
                     use_cache:bool=False,
                     work_dir:str=None
) -> Ligand:
    """Implementation function that places a specified ligand into a structure using the alphafill 
    algorithm (https://doi.org/10.1038/s41592-022-01685-y). 

    Args:
        ligand:
        simliarity_metric:
        similaritY_cutoff:
        minimize:
        min_iter:
        use_cache:
        work_dir:

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
        sim_func = ligand_mcs
    else:
        assert False
#    for i, row in df.iterrows():
#        numerator:int=len(row.atom_names.intersection(ligand_atoms))
#        denominator:int=max(len(ligand_atoms), len(row.atom_names))
#        similarity_score.append(numerator/denominator)
    df['similarity_score'] = df.apply(lambda row: sim_func( ligand, row.mol ), axis=1 )
    print(df)
    exit( 0 )
    similarity_mask = df.similarity_score >= similarity_cutoff

    if not similarity_mask.sum():
        _LOGGER.error(f"No transplants exist with similarity_cutoff of {similarity_cutoff:.3f}! Exiting...")
        exit( 1 )
    
    df = df[similarity_mask].reset_index(drop=True)
    similarity = df.similarity_score.to_numpy()
    df = df[np.isclose(similarity, np.max(similarity))].reset_index(drop=True)

    df.sort_values(by='clash_count',inplace=True)
    selected_id:str=df.asym_id[0]

    molfile:str=f"{work_dir}/temp.pdb"
    reactant:str=f"{work_dir}/temp_ligand.mol2"
    template:str=f"{work_dir}/template_ligand.mol2"
    outfile:str=f"{work_dir}/aligned_ligand.mol2"
    to_delete.extend([molfile, reactant, template, outfile, structure_start, filled_structure, transplant_info_fpath])

    parser.save_structure(molfile, stru)

    Mol2Parser().save_ligand(reactant, ligand)
    
    pp_chains:List[str] = list()
    for cc in stru.chains:
        if cc.is_polypeptide():
            pp_chains.append(f"chain {cc.name}")

    interface.pymol.general_cmd(session, [('delete', 'all'), ("load", filled_structure), ("remove", "solvent"), ("load", molfile),
                                          ("align", f"{Path(filled_structure).stem} and ({' or '.join(pp_chains)})", f"{Path(molfile).stem} and ({' or '.join(pp_chains)})"),
                                          ("delete", Path(molfile).stem),
                                          ("save", template, f"chain {selected_id} and segi {selected_id}")])
    
    analogue_code = df.analogue_id[0]
    session = interface.pymol.new_session()

    df = interface.pymol.collect(session, template, 'name x y z'.split())
    
    session = interface.pymol.new_session()
    args = [('fetch', analogue_code), ('remove', 'hydrogens')]

    for i, row in df.iterrows():
        args.append(('alter_state', 1, f"name {row['name']}", f"(x, y, z) = ({row.x}, {row.y}, {row.z})"))

    args.append(('save', template))

    interface.pymol.general_cmd(session, args)

    if ligand.is_ligand():
        outfile = mimic_torsions(template, reactant)
    else:
        fs.safe_mv(template, outfile)

    aligned_ligand:Ligand=Mol2Parser().get_ligand(outfile)
    
    for a1 in ligand.atoms:
        for a2 in aligned_ligand.atoms:
            if a1.name == a2.name:
                a1.coord = a2.coord



def _seed_ligand_analog(stru, ligand, work_dir, **kwargs):

    analog_file:str=f"{work_dir}/analog.mol2"
    analog_template_file:str=f"{work_dir}/analog_template.mol2"

    parser = Mol2Parser()
    parser.save_ligand(analog_file, kwargs['analog'])
    parser.save_ligand(analog_template_file, kwargs['analog_template'])

    analog = mimic_torsions(analog_template_file, analog_file)
    
    analog = parser.get_ligand(analog_file)

    target_location = None

    for atom in analog.atoms:
        if atom.name == kwargs['analog_atom']:
            target_location = np.array(atom.coord)
            break
    
    assert target_location is not None

    current_location = None
    for atom in ligand.atoms:
        if atom.name == kwargs['seed_atom']:
            current_location = np.array(atom.coord)
            break

    assert current_location is not None
    shift = target_location - current_location

    ligand.shift( shift )



def _seed_ligand_mole2( stru:Structure, ligand:Ligand, work_dir:str, **kwargs ) -> None:
    """
    Args:
        stru: 
        ligand:
        work_dir:

    Returns:
        Nothing.
    """
    delta = kwargs.get('delta', 0.75) 
    constraints = kwargs.get('constraints', [])
    assert constraints
    if not work_dir:
        work_dir = config['system.SCRATCH_DIR']

    relevant_constraints: List[str] = list()
    
    temp_pdb:str=f"{work_dir}/__temp.pdb"
    to_delete:List[str] = [ temp_pdb ]

    parser = PDBParser()
    parser.save_structure( temp_pdb, stru )
    cavities:List[Mole2Cavity] = interface.mole2.identify_cavities(temp_pdb)

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
        #TODO(CJ): this is where I put the actual energy/constraint evaluation
    
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

    score_mask = np.isclose(scores, np.min(scores)*2)

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
    

def minimize_ligand_only(ligand:Ligand, n_iter:int, work_dir:str) -> None:
    """Uses Rosetta to minimize only the supplied Ligand(). Designed to remove clashes,
    not perform a rigorous optimization of the molecule.

    Args:
        ligand: Ligand() to minimize.
        n_iter: How many iter's/trajectories should be used?
        work_dir: Where should work be done? Note: actually done in <work_dir>/minimize/

    Returns:
        Nothing.
    """
    sele = f'(chain {ligand.parent.name} and resi {ligand.idx})'
    geo_minimize(stru, movemap=[
        {'sele':sele, 'bb':'true', 'chi':'true', 'bondangle':'true'}, 
        {'sele':f'not ({sele})', 'bb':'false', 'chi':'false', 'bondangle':'false'}, 
    
    ],n_iter=n_iter, work_dir=f"{work_dir}/minimize/")

