"""Provides functionality for placing a ligand into a protein. All functionalty should be accessed through 
the place_ligand() function. Allows for placement of ligands via one of the below methods:

    + alphafill: 
    + rosetta_ligand: 

Author: Chris Jurich <chris.jurich@vanderbit.edu>
Date: 2023-09-14
"""

from typing import List, Tuple, Dict

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
from enzy_htp.structure import PDBParser, Structure
from enzy_htp._interface import RosettaCst

from .align_ligand import align_ligand
from .constrained_docking import constrained_docking


def place_ligand(molfile: str,
                 ligand: str,
                 code: str = None,
                 method: str = None,
                 new_res_key: Tuple = None,
                 work_dir: str = None,
                 outfile: str = None,
                 constraints: List[RosettaCst] = None,
                 **kwargs) -> Tuple[str, str]:
    """Algorithm for placing a ligand into  

    Args:
        molfile:
        ligand:
        code:
        method:
        new_res_key:
        work_dir:
        outfile:
        constraints:


    Returns:
        A Tuple[str,str] with the layout (enzyme-ligand complex in .pdb format, exact ligand conformation in .mol2 format).
    """
    #TODO(CJ): check for None values and throw errors if they are given
    if not work_dir:
        work_dir = config['system.SCRATCH_DIR']

    fs.safe_mkdir(work_dir)

    placed_ligand: str = None

    if method == "alphafill":
        placed_ligand = _place_alphafill(molfile, ligand, code=code, work_dir=work_dir, **kwargs)
    elif method == "rosetta_ligand":
        place_ligand = _place_rosetta_ligand(molfile, ligand, new_res_key=new_res_key, work_dir=work_dir, constraints=constraints, **kwargs)
        pass
    else:
        _LOGGER.error(
            f"The supplied method placement method {method} is not supported. Allowed methods are 'alphafill' and 'rosetta_ligand'. Exiting..."
        )
        exit(1)

    session = interface.pymol.new_session()

    args = [('load', placed_ligand)]

    if new_res_key is not None:
        args.extend([
            ('alter', 'all', f'chain="{new_res_key[0]}"'),
            ('alter', 'all', f'resi="{new_res_key[1]}"'),
            ('alter', 'all', f'resn="{new_res_key[2]}"'),
            ('alter', 'all', 'segi=""'),
        ])

    if not outfile:
        temp_path = Path(molfile)
        outfile = f"{work_dir}/{temp_path.stem}_placed{temp_path.suffix}"

    args.extend([('save', placed_ligand), ('load', molfile), ('save', outfile)])

    interface.pymol.general_cmd(session, args)

    return (outfile, placed_ligand)


def _place_rosetta_ligand(molfile: str, reactant: str, new_res_key: Tuple, work_dir: str, constraints: List[RosettaCst], **kwargs) -> str:
    """TODO(CJ)

    Args:

    Returns:
        
    """
    if not work_dir:
        work_dir = config['system.SCRATCH_DIR']

    relevant_constraints: List[str] = list()

    for cc in constraints:
        if cc.contains(new_res_key[0], new_res_key[1]):
            relevant_constraints.append(cc)

    if not relevant_constraints:
        _LOGGER.error(
            f"To place a ligand with RosettaLigand, you must supply at least one constraint involving the target ligand. Exiting...")
        exit(1)
    else:
        _LOGGER.info(f"Found {len(relevant_constraints)} constraints that involve the ligand! Continuing...")

    temppath = Path(molfile)
    other = str(temppath.parent / "other.pdb")

    session = interface.pymol.new_session()
    interface.pymol.general_cmd(
        session,
        [
            ('load', reactant),
            ('alter', 'all', f"chain='{new_res_key[0]}'"),
            ('alter', 'all', f"resi='{new_res_key[1]}'"),
            ('alter', 'all', f"resn='{new_res_key[2]}'"),
            ('alter', 'all', f"segi=''"),
            #('load', molfile),
            ('save', reactant),
            ('delete', 'all')
        ])
    interface.pymol.general_cmd(session, [('load', molfile)])

    df: pd.DataFrame = interface.pymol.collect(session, 'memory', 'chain resn resi'.split())
    need_params: List[str] = list()

    for i, row in df.iterrows():
        key = (row.chain, row.resn, row.resi)

        if key in need_params:
            continue

        if key[1].upper() in chem.METAL_CENTER_MAP:
            continue

        if key[1].upper() in chem.THREE_LETTER_AA_MAPPER:
            continue

        need_params.append(key)

    param_files: List[str] = list()
    for needs in need_params:
        temp_file: str = f"{work_dir}/temp.mol2"
        interface.pymol.general_cmd(session, [("save", temp_file, f"chain {needs[0]} and resn {needs[1]} and resi {needs[2]}")])
        param_files.append(interface.rosetta.parameterize_ligand(temp_file, needs[1])[0])

    param_files.append(interface.rosetta.parameterize_ligand(reactant, new_res_key[2])[0])
    parser = eh.PDBParser()

    interface.pymol.general_cmd(session, [('load', reactant), ('save', other)])

    structure = parser.get_structure(other)
    #TODO(CJ): make this less error prone
    #TODO(CJ): need to look through and find all the parts where there are
    constrained_docking(structure, [reactant], param_files, relevant_constraints, work_dir=work_dir)
    #what needs to happen here?
    # 1. figure out which constraints involve the reactant we care about
    # 2. create the input files that we'll need
    # 3. run the rosetaLigand run
    # 4. select placement with best constraint score
    assert False


def _create_xml(work_dir: str, reactants: List[str], chain_names: List[str], use_cache: bool) -> str:
    """

    Args:
        

    Returns:
        The relative filepath of the RosettaScripts .xml file.
    """
    reactant_dicts = list()
    for ridx, rr in enumerate(reactants):
        temp = dict()
        session = interface.pymol.new_session()
        df: pd.DataFrame = interface.pymol.collect(session, rr, "chain resn x y z name".split(), sele='not elem H')
        temp['resn'] = df.resn[0]
        temp['chain'] = chain_names[ridx]
        temp['volume'] = interface.rdkit.volume(rr)
        points = np.transpose(np.array([df.x.to_numpy(), df.y.to_numpy(), df.z.to_numpy()]))
        distances = list()
        n_points = len(points)
        for p1 in range(n_points):
            for p2 in range(n_points):
                distances.append(np.sqrt(np.sum((points[p1] - points[p2])**2)))
        temp['length'] = np.max(np.array(distances))
        reactant_dicts.append(temp)

    reactant_dicts.sort(key=lambda dd: dd['chain'])

    fname: str = f"{work_dir}/__script.xml"
    if not use_cache:
        fs.safe_rm(fname)

    fpath = Path(fname)
    if fpath.exists() and use_cache:
        _LOGGER.info(f"Using cached RosettaScripts .xml file at {fpath.absolute()}")
        return fname

    elements: List[Dict] = [
        {
            'parent': 'SCOREFXNS',
            'tag': 'ScoreFunction',
            'name': 'ligand_soft_rep',
            'weights': 'ligand_soft_rep'
        },
        {
            'parent': 'SCOREFXNS',
            'tag': 'ScoreFunction',
            'name': 'hard_rep',
            'weights': 'ligand'
        },
        {
            'parent': 'SCOREFXNS.ScoreFunction',
            'tag': 'Reweight',
            'scoretype': 'coordinate_constraint',
            'weight': '1.0'
        },
        {
            'parent': 'SCOREFXNS.ScoreFunction',
            'tag': 'Reweight',
            'scoretype': 'atom_pair_constraint',
            'weight': '1.0'
        },
        {
            'parent': 'SCOREFXNS.ScoreFunction',
            'tag': 'Reweight',
            'scoretype': 'angle_constraint',
            'weight': '1.0'
        },
        {
            'parent': 'SCOREFXNS.ScoreFunction',
            'tag': 'Reweight',
            'scoretype': 'dihedral_constraint',
            'weight': '1.0'
        },
        {
            'parent': 'SCOREFXNS.ScoreFunction',
            'tag': 'Reweight',
            'scoretype': 'chainbreak',
            'weight': '1.0'
        },
    ]

    for rd in reactant_dicts:
        cc = rd['chain']
        elements.extend([
            {
                'parent': 'LIGAND_AREAS',
                'tag': 'LigandArea',
                'name': f'docking_sidechain_{cc.lower()}',
                'chain': f'{cc.upper()}',
                'cutoff': '6.0',
                'add_nbr_radius': 'true',
                'all_atom_mode': 'false'
            },
            {
                'parent': 'LIGAND_AREAS',
                'tag': 'LigandArea',
                'name': f'final_sidechain_{cc.lower()}',
                'chain': f'{cc.upper()}',
                'cutoff': '6.0',
                'add_nbr_radius': 'true',
                'all_atom_mode': 'false'
            },
            {
                'parent': 'LIGAND_AREAS',
                'tag': 'LigandArea',
                'name': f'final_backbone_{cc.lower()}',
                'chain': f'{cc.upper()}',
                'cutoff': '7.0',
                'add_nbr_radius': 'false',
                'all_atom_mode': 'true',
                'Calpha_restraints': "0.3"
            },
            {
                'parent': 'INTERFACE_BUILDERS',
                'tag': 'InterfaceBuilder',
                'name': f'side_chain_for_docking_{cc.lower()}',
                'ligand_areas': f'docking_sidechain_{cc.lower()}'
            },
            {
                'parent': 'INTERFACE_BUILDERS',
                'tag': 'InterfaceBuilder',
                'name': f'side_chain_for_final_{cc.lower()}',
                'ligand_areas': f'final_sidechain_{cc.lower()}'
            },
            {
                'parent': 'INTERFACE_BUILDERS',
                'tag': 'InterfaceBuilder',
                'name': f'backbone_{cc.lower()}',
                'ligand_areas': f'final_backbone_{cc.lower()}',
                'extension_window': '3'
            },
            {
                'parent': 'MOVEMAP_BUILDERS',
                'tag': 'MoveMapBuilder',
                'name': f'docking_{cc.lower()}',
                'sc_interface': f'side_chain_for_docking_{cc.lower()}',
                'minimize_water': 'false'
            },
            {
                'parent': 'MOVEMAP_BUILDERS',
                'tag': 'MoveMapBuilder',
                'name': f'final_{cc.lower()}',
                'sc_interface': f'side_chain_for_final_{cc.lower()}',
                'bb_interface': f'backbone_{cc.lower()}',
                'minimize_water': 'false'
            },
        ])

    for ridx, rd in enumerate(reactant_dicts):
        cc = rd['chain']
        rd['grid_name'] = f"grid_{ridx+1}"

        elements.extend([
            {
                'parent': 'ROSETTASCRIPTS',
                'tag': 'SCORINGGRIDS',
                'ligand_chain': cc.upper(),
                'width': str(rd['length'] * 2),
                'append_elements_only': True,
                'name': rd['grid_name'],
                'child_nodes': [
                    {
                        'parent': 'SCORINGGRIDS',
                        'tag': 'ClassicGrid',
                        'grid_name': 'classic',
                        'weight': '1.0'
                    },
                    {
                        'parent': 'SCORINGGRIDS',
                        'tag': 'HbdGrid',
                        'grid_name': 'hbd_grid',
                        'weight': '1.0'
                    },
                    {
                        'parent': 'SCORINGGRIDS',
                        'tag': 'HbaGrid',
                        'grid_name': 'hba_grid',
                        'weight': '1.0'
                    },
                ]
            },
        ])

    elements.extend([{'parent': 'MOVERS', 'tag': 'AddOrRemoveMatchCsts', 'name': 'cstadd', 'cst_instruction': 'add_new'}])

    for rd in reactant_dicts:
        cc: str = rd['chain']
        radius: float = rd['length']
        elements.extend([{
            'parent': 'MOVERS',
            'tag': 'Transform',
            'name': f'transform_{cc.lower()}',
            'chain': f'{cc.upper()}',
            'box_size': str(1.5 * radius),
            'move_distance': '0.1',
            'angle': '45',
            'cycles': '200',
            'repeats': '10',
            'temperature': '25',
            'grid_set': rd['grid_name']
        }, {
            'parent': 'MOVERS',
            'tag': 'HighResDocker',
            'name': f'high_res_docker_{cc.lower()}',
            'cycles': '6',
            'repack_every_Nth': '3',
            'scorefxn': 'ligand_soft_rep',
            'movemap_builder': f'docking_{cc.lower()}'
        }, {
            'parent': 'MOVERS',
            'tag': 'FinalMinimizer',
            'name': f'final_{cc.lower()}',
            'scorefxn': 'hard_rep',
            'movemap_builder': f'final_{cc.lower()}'
        }])

    elements.append({
        'parent': 'MOVERS',
        'tag': 'InterfaceScoreCalculator',
        'name': 'add_scores',
        'chains': ','.join(map(lambda ss: ss['chain'].upper(), reactant_dicts)),
        'scorefxn': 'hard_rep'
    })

    for rd in reactant_dicts:
        cc: str = rd['chain']

        elements.extend([{
            'parent': 'MOVERS',
            'tag': 'ParsedProtocol',
            'name': f'low_res_dock_{cc.lower()}',
            'child_nodes': [
                {
                    'tag': 'Add',
                    'mover_name': 'cstadd'
                },
                {
                    'tag': 'Add',
                    'mover_name': f'transform_{cc.lower()}'
                },
            ]
        }])

    elements.extend([
        {
            'parent': 'MOVERS',
            'tag': 'ParsedProtocol',
            'name': 'high_res_dock',
            'child_nodes': [deepcopy({
                'tag': 'Add',
                'mover_name': f'high_res_docker_{rd["chain"].lower()}'
            }) for rd in reactant_dicts]
        },
        {
            'parent': 'MOVERS',
            'tag': 'ParsedProtocol',
            'name': 'reporting',
            'child_nodes': [{
                'tag': 'Add',
                'mover_name': 'add_scores'
            }]
        },
    ])

    elements.extend(
        [deepcopy({
            'parent': 'PROTOCOLS',
            'tag': 'Add',
            'mover_name': f'low_res_dock_{rd["chain"].lower()}'
        }) for rd in reactant_dicts] + [{
            'parent': 'PROTOCOLS',
            'tag': 'Add',
            'mover_name': 'high_res_dock'
        }, {
            'parent': 'PROTOCOLS',
            'tag': 'Add',
            'mover_name': 'reporting'
        }])

    interface.rosetta.write_script(fname, elements)

    _LOGGER.info(f"Saved new RosettaScripts .xml file at {fpath.absolute()}!")
    return fname


def _place_alphafill(molfile: str,
                     reactant: str,
                     code: str,
                     work_dir: str,
                     similarity_cutoff: float = 0.75,
                     clash_radius: float = 2.0,
                     **kwargs) -> str:
    """Implementation function that places a specified ligand into a molfile structure using the alphafill 
    algorithm (https://doi.org/10.1038/s41592-022-01685-y). Assumes that the input ligand has PDB style naming 
    and is in the .mol2 format. When scanning alphafill transplants, the algorithm looks for exact atom label
    matches. If no exact matches exist, then it looks for candidates where the reactant's atom names are a subset 
    of the candidates. If no satisfactory matches are found, it looks for candidates where the candidate's atom
    names are a subset of the reactants. For each subcategory, candidates template molecules are ranked by the 
    number of clashes with the protein structure, with fewer being better. The satisfactory template molecule is
    saved to a .mol2 file in the work_dir and returned. If no template is found, an empty str() is returned.


    Args:
        molfile:
        reactant:
        code:
        work_dir:
        similarity_cutoff:
        clash_radius:

    Returns:
        THe path to the placed reactant in a .mol2 file.

    """
    #TODO(CJ): get the kwargs unpacked here
    _LOGGER.info(f"Beginnning placement of ligand {reactant} into {molfile}...")
    _LOGGER.info(f"Calling out to AlphaFill to fill structure...")
    filled_structure: str = interface.pymol.fetch(code, work_dir)
    filled_structure = interface.alphafill.fill_structure(filled_structure)
    _LOGGER.info("Filled structure using AlphaFill!")

    session = interface.pymol.new_session()
    interface.pymol.general_cmd(session, [("load", filled_structure), ("remove", "solvent"), ("load", molfile),
                                          ("align", Path(filled_structure).stem, Path(molfile).stem)])
    f_df: pd.DataFrame = interface.pymol.collect(session,
                                                 'memory',
                                                 "chain resn resi segi name elem x y z".split(),
                                                 sele=f"{Path(filled_structure).stem} and (not elem H)")
    p_df: pd.DataFrame = interface.pymol.collect(session,
                                                 'memory',
                                                 "chain resn resi segi name elem x y z".split(),
                                                 sele=f"{Path(molfile).stem} and (not elem H)")
    interface.pymol.general_cmd(session, [("delete", Path(molfile).stem)])

    protein_tks = list()
    candidates = list()
    for (chain, resi, resn, segi) in set(zip(f_df.chain, f_df.resi, f_df.resn, f_df.segi)):
        if resn in eh.chemical.THREE_LETTER_AA_MAPPER:
            protein_tks.append(f"( chain {chain} and resi {resi} and resn {resn} and segi {segi})")
            continue

        if resn.upper() in "MG ZN CL NA".split():  #TODO(CJ): make this a function in chemical
            continue
        temp = f_df[(f_df.chain == chain) & (f_df.resi == resi) & (f_df.resn == resn) & (f_df.segi == segi)].reset_index(drop=True)
        candidates.append(deepcopy({'df': temp, 'names': set(temp.name.to_list()), 'key': (chain, resi, resn, segi)}))

    protein_sele: str = "not elem H and " + " or ".join(protein_tks)
    interface.pymol.general_cmd(session, [("delete", "all"), ("load", filled_structure), ("remove", "solvent")])
    #p_df:pd.DataFrame=interface.pymol.collect(session, 'memory', "chain resn resi segi name elem x y z".split(), sele=protein_sele)

    for cc in candidates:
        cc['clashes'] = _count_clashes(p_df, cc['df'], clash_radius)

    interface.pymol.general_cmd(session, [("delete", "all"), ("load", reactant), ("remove", "solvent")])
    r_df: pd.DataFrame = interface.pymol.collect(session, 'memory', "chain resn resi segi name elem x y z".split(), sele="not elem H")

    rct = {'df': r_df, 'names': set(r_df.name.to_list())}

    exact_matches = list()
    rct_is_subset = list()
    can_is_subset = list()

    for cc in candidates:
        if rct['names'] == cc['names']:
            exact_matches.append(cc)
            continue

        if rct['names'].issubset(cc['names']):
            rct_is_subset.append(cc)
            continue

        if cc['names'].issubset(rct['names']):
            can_is_subset.append(cc)
            continue

    template = str()
    if exact_matches:
        exact_matches.sort(key=lambda dd: dd['clashes'])
        row = exact_matches[0]['df'].iloc[0]

        _LOGGER.info(f"Found exact match at {row.chain}.{row.resi}.{row.resn}! Transplanting")
        temp_path = Path(reactant)
        outfile = temp_path.parent / f"{temp_path.stem}_placed.mol2"
        template: str = f"{work_dir}/template.mol2"
        interface.pymol.general_cmd(
            session, [("delete", "all"), ("load", filled_structure), ("load", molfile),
                      ("align", Path(filled_structure).stem, Path(molfile).stem),
                      ("save", template,
                       f"{Path(filled_structure).stem} and chain {row.chain} and resn {row.resn} and resi {row.resi} and segi {row.segi}")])
        return align_ligand(template, reactant, outfile=outfile)

    assert False
    if rct_is_subset:
        pass


def _count_clashes(df1: pd.DataFrame, df2: pd.DataFrame, cutoff: float) -> int:
    """Counts the number of clashes between the points in df1 and df2

    Args:
        df1:
        df2:
        cutoff:

    Returns:
        pass
    """
    #TODO(CJ): do checks on the input pandas dataframes
    points1 = np.transpose(np.array([df1.x.to_numpy(), df1.y.to_numpy(), df1.z.to_numpy()]))
    points2 = np.transpose(np.array([df2.x.to_numpy(), df2.y.to_numpy(), df2.z.to_numpy()]))

    count: int = 0

    for pp1 in points1:
        for pp2 in points2:
            count += np.sqrt(np.sum((pp1 - pp2)**2)) <= cutoff

    return count
