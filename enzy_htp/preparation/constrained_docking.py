"""Drive code for the constrained docking protocol in enzy_htp. This implementation makes use of RosettaLigand and the enzdes constraints
to dock structures with constrained geometries. Users are meant to call the `constrained_docking()` method and none of the other functions
found in this sub-module. Consult `constrained_docking()` for more information on the protocol and its settings. This function is used in 
both enzy_htp/preparation/reactive_docking.py and enzy_htp/preparation/place_ligand.py with differing settings.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-09-20
"""
import os
from pathlib import Path
from typing import List, Tuple
from copy import deepcopy

import pandas as pd
import numpy as np

from enzy_htp.core import file_system as fs
from enzy_htp.core import _LOGGER
from enzy_htp._interface.rosetta_interface import RosettaCst
from enzy_htp.structure import PDBParser, Structure

from enzy_htp import interface, config


def constrained_docking(structure: Structure,
                        reactants: List[str],
                        param_files: List[str],
                        constraints: List[RosettaCst],
                        work_dir: str = None,
                        use_cache: bool = False,
                        rng_seed: int = 1996,
                        n_struct: int = 100,
                        minimize: bool = True) -> pd.DataFrame:
    """Performs constrained docking using the RosettaLigand docking algorithm and the enzdes constraints system available within Rosetta.
    All results are contained within a pandas DataFrame format returned by the function. The DataFrame contains what would normally be 
    found in the score.sc file.
    
    Args:
        structure: The Structure object serving as the basis for the docking.
        reactants: A List[str] specifying paths of the reactants to perform docking on. 
        param_files: A List[str] of .params files for ligands in the system.
        constraints: A List[RosettaCst] defining geometric constraints for the system.
        work_dir: Where all temp files should be saved. Optional, defaults to config['system.SCRATCH_DIR'].
        use_cache: Should the cache be used? Defaults to False.
        rng_seed: int() rng seed to use. Defaults to 1996.
        n_struct: Number of docked conformers to produce. Defaults to 100.
        minimize: Should side chain minimization be performed? Defaults to True.
        
    Returns:
        A pandas Dataframe containing all results for the docking run.         
    """
    if work_dir is None:
        work_dir = config['system.SCRATCH_DIR']


    (start_pdb, cst_file) = _integrate_csts(structure, constraints, work_dir, use_cache)

    xml_file: str = _create_xml(work_dir, reactants, use_cache, minimize) 

    options_file: str = _make_options_file(start_pdb, xml_file, param_files, work_dir, rng_seed, n_struct, use_cache, cst_file)

    df: pd.DataFrame = _docking_run(options_file, use_cache)

    return df


def _integrate_csts(stru: Structure, csts: List[RosettaCst], work_dir: str, use_cache: bool) -> Tuple[str, str]:
    """Given a structure and constraints, combines all necessary information so that an enzdes-constrained RosettaLigand run
    can be performed. Results in a .pdb file with the appropriate REMARK lines and a .cst constraints file that documents
    the described constraints.
    
    Args:
        stru: The Structure to save into a .pdb file and apply constraints to.  
        csts: A list() of RosettaCst's to incorporate into the .pdb and .cst files. 
        work_dir: The directory to save the .pdb and .cst files in. 
        use_cache: Should existing files be used when available? 

    Returns:
        A Tuple() with the format (pdb file, .cst file).

    """
    _LOGGER.info("Beginning RosettaCst constraint integration...")
    parser = PDBParser()
    file_str = parser.get_file_str(stru, if_renumber=False, if_fix_atomname=False)

    pdb_content: List[str] = ["HEADER                                            xx-MMM-xx"]
    cst_content: List[str] = list()
    for cidx, cst in enumerate(csts):
        pdb_content.append(cst.create_pdb_line(cidx + 1))
        cst_content.extend(cst.create_cst_lines())

    pdb_file: str = f"{work_dir}/start.pdb"
    cst_file: str = f"{work_dir}/rdock.cst"

    if not use_cache:
        fs.safe_rm(pdb_file)
        fs.safe_rm(cst_file)

    if not Path(pdb_file).exists():
        fs.write_lines(pdb_file, pdb_content + file_str.splitlines())

    if not Path(cst_file).exists():
        fs.write_lines(cst_file, cst_content)

    _LOGGER.info("RosettaCst constraint integration successful! Relevant files:")
    _LOGGER.info(f"\t.pdb file: {Path(pdb_file).absolute()}")
    _LOGGER.info(f"\t.cst file: {Path(cst_file).absolute()}")

    return (pdb_file, cst_file)


def _create_xml(work_dir: str, reactants: List[str], use_cache: bool, minimize:bool) -> str:
    """Creates the RosettaScripts .xml file for the RosettaLigand run. Can be run in the cache or minmization mode 
    which uses the existing .xml script and allow for side chain-ligand minimization respectively , respectively. 
    This function assumes that the reactants are in a format where their associated chain id and residue names are correct.
    

    Args:
        work_dir: The dir as a str() where the .xml script should be saved.
        reactants: A List[str] of reactant names. 
        use_cache: Should existing files be used when possible?
        minimize: Should minimization be performed.
        

    Returns:
        The relative filepath of the RosettaScripts .xml file.
    """
    fname: str = f"{work_dir}/__script.xml"
    if not use_cache:
        fs.safe_rm(fname)
    elif fs.has_content(fname):
        _LOGGER.info(f"Using cached RosettaScripts .xml file at {fname}")
        return fname



    reactant_dicts = list()
    for rr in reactants:
        temp = dict()
        session = interface.pymol.new_session()
        df: pd.DataFrame = interface.pymol.collect(session, rr, "chain resn x y z name".split(), sele='not elem H')
        temp['resn'] = df.resn[0]
        temp['chain'] = df.chain[0]
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

    _LOGGER.info(f"Saved new RosettaScripts .xml file at {fname}!")
    return fname


def _make_options_file(pdb_file: str,
                      xml_file: str,
                      param_files: List[str],
                      work_dir: str,
                      rng_seed: int,
                      n_struct: int,
                      use_cache: bool,
                      cst_file: str = None) -> str:
    """Makes the options.txt file that the docking run will actually use. Takes a variety of arguments and 
    can used cached values if needed. This function DOES NOT make any checks to the inputs.
    
    Args:
        pdb_file: The .pdb file (with constraints) to use. 
        xml_file: The validated RosettaScripts .xml file to be used for docking.
        param_files: The list() of reactant .params files.
        work_dir: The working directory 
        rng_seed: rng seed to be used during dcking.
        n_struct: Number of strutures to make as an int().
        use_cache: Should we used cached values? 
        cst_file: The contraints file to be used. Optional. 

    Returns:
        Path to the options.txt file with all the Rosetta options.

    """
    _LOGGER.info("Beginning creation of options file for Rosetta...")
    content: List[str] = [
        "-keep_input_protonation_state",
        "-run:constant_seed",
        f"-run:jran {int(rng_seed)}",
        "-in:file",
        f"    -s '{Path(pdb_file).name}'",
    ]

    for pf in param_files:
        content.append(f"    -extra_res_fa '{Path(pf).absolute()}'")

    stub_parent: str = os.path.expandvars(
        f"${config['rosetta.ROSETTA3']}/database/chemical/residue_type_sets/fa_standard/residue_types/protonation_states/")
    for stub in "GLU_P1.params GLU_P2.params LYS_D.params ASP_P1.params TYR_D.params HIS_P.params ASP_P2.params".split():
        content.append(f"    -extra_res_fa '{stub_parent}/{stub}'")

    content.extend([
        "-run:preserve_header",
        "-packing",
        "    -ex1",
        "    -ex2aro",
        "    -ex2 ",
        "    -no_optH false",
        "    -flip_HNQ true",
        "    -ignore_ligand_chi true",
    ])

    if cst_file:
        content.extend([f"-enzdes:cstfile '{Path(cst_file).absolute()}'"])
    else:
        _LOGGER.warning("No constraints supplied! This will hurt reaction complex accuracy!")

    content.extend([
        "-parser",
        f"   -protocol {Path(xml_file).absolute()}",
        "-out",
        f"   -file:scorefile 'score.sc'",
        "   -level 200",
        f"   -nstruct {n_struct}",
        "   -overwrite",
        "   -path",
        f"       -all './complexes'",
    ])

    qsar_grid: str = str(Path(f"{work_dir}/complexes/qsar_grids/").absolute())
    content.append(f"-qsar:grid_dir {qsar_grid}")

    fname = Path(work_dir) / "options.txt"
    score_file: str = f"{work_dir}/complexes/score.sc"

    if not use_cache:
        _LOGGER.info("Not using cache mode. Deleting the following (if they already exist):")
        _LOGGER.info(f"\toptions file: {fname}")
        _LOGGER.info(f"\tscore file: {score_file}")
        _LOGGER.info(f"\tenzyme-reactant complexes directory: {work_dir}/complexes")
        _LOGGER.info(f"\tqsar_gird directory: {qsar_grid}")
        fs.safe_rm(fname)
        fs.safe_rm(score_file)
        fs.safe_rmdir(f"{work_dir}/complexes/")
        fs.safe_rmdir(qsar_grid)

    fs.safe_mkdir(f"{work_dir}/complexes/")
    fs.safe_mkdir(qsar_grid)

    if not fname.exists():
        _LOGGER.info(f"Wrote the below settings to {fname}:")
        for ll in content:
            _LOGGER.info(f"\t{ll}")
        fs.write_lines(fname, content)
    else:
        _LOGGER.info(f"Cache mode enabled. Using below settings from {fname}:")
        content: List[str] = fs.lines_from_file(fname)
        for ll in content:
            _LOGGER.info(f"\t{ll}")

    option_file = fname.absolute()

    return str(fname)


def _docking_run(option_file: str, use_cache: bool) -> pd.DataFrame:
    """Executes the actual docking run prepared by the rest of the functions in this sub-module. Changes directories
    to where the options.txt file is located. All results are returned in a pandas DataFrame.

    Args:
        option_file: Name of the options.txt file being used as a str().
        use_cache: Should we use existing results when available.

    Returns:
        The results in a pandas DataFrame format.
    """
    _LOGGER.info("Beginning RosettaLigand docking run...")
    opt_path = Path(option_file)

    if use_cache:
        csv_file = Path(str(opt_path.parent / "scores.csv"))
        _LOGGER.info(f"Cache mode enabled. Looking for {csv_file} ...")
        if csv_file.exists():
            df = pd.read_csv(csv_file)
            df['selected'] = True
            _LOGGER.info(f"Found file {csv_file}. Checking for existence of {len(df)} output .pdbs")

            for i, row in df.iterrows():
                if (i + 1) % 25 == 0:
                    _LOGGER.info(f"\tChecking {i+1} of {len(df)}...")
                if not Path(row.description).exists():
                    _LOGGER.info(f"\t{row.description} is missing. Cannot use cached results!")
                    break
            else:
                _LOGGER.info("All cached .pdbs exist! Using cached RosettaLigand structures!")
                return df
        else:
            _LOGGER.info(f"{csv_file} not found. Continuing with standard run")

    _LOGGER.info("Beginning RosettaLigand geometry sampling step...")
    start_dir: str = os.getcwd()
    os.chdir(str(opt_path.parent))

    interface.rosetta.run_rosetta_scripts([f"@{opt_path.name}"])

    os.chdir(start_dir)

    df: pd.DataFrame = interface.rosetta.parse_score_file(str(opt_path.parent / "complexes/score.sc"))

    df['description'] = df.apply(lambda row: f"{opt_path.parent}/complexes/{row.description}.pdb", axis=1)

    df['selected'] = True

    df.to_csv(str(opt_path.parent / "scores.csv"), index=False)

    _LOGGER.info("Completed RosettaLigand geometry sampling!")

    return df
