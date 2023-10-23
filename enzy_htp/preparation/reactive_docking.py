"""Driver for the reactive docking functionality available in enzy_htp. The only function that that should be 
called is dock_reactants(). All others are implementation functions that SHOULD NOT be used. 

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-07-28
"""
import os
from pathlib import Path
from collections import defaultdict
from typing import List, Tuple, Set, Dict
from copy import deepcopy

import numpy as np
import pandas as pd

from enzy_htp import interface, config, _LOGGER, binding_energy
from enzy_htp._interface.rosetta_interface import RosettaCst
import enzy_htp.structure.structure_operation as stru_oper

import enzy_htp.chemical as chem
from enzy_htp import mutation as mm
from enzy_htp.structure import PDBParser, Mol2Parser, Structure, Ligand
from enzy_htp.core import file_system as fs

#TODO(CJ): I should go through and use existing values if cache mode is selected
#TODO(CJ): need to be able to specify the size of the grid

def dock_reactants(structure: Structure,
                   constraints: List[RosettaCst] = None,
                   n_struct: int = 1000,
                   cst_energy: float = None,
                   clash_distance: float = 2.0,
                   clash_cutoff: int = 1,
                   max_sasa_ratio:float = 0.60,
                   cluster_distance: float = 2.0,
                   use_cache: bool = True,
                   rng_seed: int = 1996,
                   work_dir: str = None,
                   save_work_dir: bool = True
                   ) -> Structure:
    """Given a structure containing reactants and constraints defining a pre-reaction complex, perform reactive docking
    to produce the desired geometry. Function consists of the following steps:

    1. system validation
    2. low-resolution geometry sampling and minimization
    3. geometry pruning
        a. SASA ratio filtering
        b. clash count filtering
        c. constraint energy filtering
        d. QM energy ranking
    4. QM-active site optimization

    Args:
        structure:
        constraints:
        n_struct:
        cst_energy:
        clash_distance:
        clash_cutoff:
        max_sasa_ratio:
        cluster_distance:
        use_cache:
        rng_seed:
        work_dir:
        save_work_dir:

    Returns:
        A Structure() with the specified, constrained geometry. 

    """

    if work_dir is None:
        work_dir = config["system.SCRATCH_DIR"]

    if cst_energy is None:
        cst_energy = len(constraints)*2000.0

    fs.safe_mkdir(work_dir)

    _log_settings(constraints, n_struct, cst_energy, clash_distance, clash_cutoff, max_sasa_ratio,
                cluster_distance, use_cache, rng_seed, work_dir, save_work_dir)
    
    _validate_system( structure, constraints )
   
    (param_files, charge_mapper) = _parameterize_system(structure, work_dir)

    df = _dock_system(structure, constraints, n_struct, param_files, charge_mapper, work_dir, rng_seed, use_cache )
    
    _evaluate_SASA(df, max_sasa_ratio)

    _evaluate_clashes(df, clash_distance, clash_cutoff)

    _evaluate_csts(df, constraints, cst_cutoff)
   
    _evaluate_qm(df, structure, charge_mapper, cluster_distance)

    final_geometry:str = _qm_minimization(df, charge_mapper, cluster_distance, constraints)
    
    parser = PDBParser()
    final_stru = parser.get_structure(final_geometry)
    stru_oper.update_residues(structure, final_stru)
   
    if not save_work_dir:
        _LOGGER.info(f"save_work_dir set to False! Deleting {work_dir}")
        fs.safe_rmdir( work_dir )

    return structure 

def _validate_system( stru:Structure, constraints:List[RosettaCst] ) -> None:
    """Checks if the supplied list of contraints exists and is compatible with the supplied Structure() object.
    If it is not, a detailed error message is logged and the program exits.
    
    Args:
        stru: The Structure() to validate the proposed constraints in.
        constraints: A List[RosettaCst] objects 

    Returns:
        Nothing.
    """
    _LOGGER.info("Validating the supplied system...")
    msg:str=''
    error:bool = False
    for cst in constraints:
        r1_found, r2_found = False, False
        for rr in stru.residues:
            if rr.chain.name == cst.rchain_1 and rr.name == cst.rname_1 and rr.idx == cst.rnum_1:
                atom_names:List[str] = [aa.name for aa in rr.atoms]
                for aa in cst.ratoms_1:
                    if aa not in atom_names:
                        break
                else:
                    r1_found = True

            if rr.chain.name == cst.rchain_2 and rr.name == cst.rname_2 and rr.idx == cst.rnum_2:
                atom_names:List[str] = [aa.name for aa in rr.atoms]
                for aa in cst.ratoms_2:
                    if aa not in atom_names:
                        break
                else:
                    r2_found = True

        if not r1_found:
            msg += f"({cst.rchain_1}.{cst.rnum_1}.{cst.rname_1}),"
            error = True


        if not r2_found:
            msg += f"({cst.rchain_2}.{cst.rnum_2}.{cst.rname_2}),"
            error = True


    if error:
        _LOGGER.error(f"Errors while trying to validate the supplied system. Could not locate residues or had incorrect atom names: {msg[:-1]}")
        _LOGGER.error("Exiting...")
        exit( 1 )
    else:   
        _LOGGER.info("Validation complete! Supplied constraints can be applied to the supplied Structure()") 

        

def _parameterize_reactant(reactant: str, reactant_name: str, conformers: str = None) -> Tuple[str, int]:
    """Parameterizes reactant for RosettaLigand docking. Adds both a conformers and formal charge
    lines to the file.
    
    Args:
        reactant: Path to the reactant.
        reactant_name: 3 letter PDB style reactant name as a str().
        conformers: Path to files containing conformers. Optional.
    
    Returns:
        Path to the Rosetta .params file for the reactant.

    """

    if len(reactant_name) != 3:
        _LOGGER.error(f"Invalid reactant name supplied: '{reactant_name}'. Must be alphanumeric and 3 letters long. Exiting...")
        exit(1)

    #TODO(CJ): create the conformers if there aren't any
    session = interface.pymol.new_session()
    (param_file, pdb_file) = interface.rosetta.parameterize_ligand(reactant, reactant_name, conformers=conformers)
    fs.safe_rm(pdb_file)

    sdf_file: str = interface.pymol.convert(session, reactant, new_ext=".sdf")
    fcharge: int = interface.bcl.calculate_formal_charge(sdf_file)

    fs.safe_rm(sdf_file)

    lines: List[str] = fs.lines_from_file(param_file)

    lines.append(f"NET_FORMAL_CHARGE {int(fcharge)}")

    fs.safe_rm(param_file)
    fs.write_lines(param_file, lines)

    return (param_file, fcharge)


def _integrate_csts(stru: Structure, csts: List[RosettaCst], work_dir: str, use_cache: bool) -> Tuple[str, str]:
    """
    
    Args:
        stru:
        csts:
        work_dir:

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


def _create_xml(stru:Structure, work_dir: str, use_cache: bool) -> str:
    """

    Args:
        

    Returns:
        The relative filepath of the RosettaScripts .xml file.
    """
    parser = Mol2Parser()
    temp_mol2:str=f"{work_dir}/__temp_mol2.mol2"
    reactant_dicts = list()
    for rr in stru.residues:
        if not rr.is_ligand():
            continue
        temp = dict()
        temp['resn'] = rr.name
        temp['chain'] = rr.chain.name
        parser.save_ligand(temp_mol2, rr )
        temp['volume'] = interface.rdkit.volume(temp_mol2)
        reactant_dicts.append( temp )
        coords = [np.array(aa.coord) for aa in rr.atoms]
        max_dist:float = -1.0
        for c1 in coords:
            for c2 in coords:
                dist = np.sqrt(np.sum((c1-c2)**2))
                if dist > max_dist:
                    max_dist = dist
        temp['length'] = max_dist
    fs.safe_rm( temp_mol2 )

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


    for ridx, rd in enumerate(reactant_dicts):
        cc = rd['chain']
        rd['grid_name'] = f"grid_{ridx+1}"

        elements.extend([
            {
                'parent': 'ROSETTASCRIPTS',
                'tag': 'SCORINGGRIDS',
                'ligand_chain': cc.upper(),
                'width': '20', 
                'append_elements_only': True,
                'name': rd['grid_name'],
                'child_nodes': [
                    {
                        'parent': 'SCORINGGRIDS',
                        'tag': 'ClassicGrid',
                        'grid_name': 'classic',
                        'weight': '1.0'
                    },
                ]
            },
        ])

    elements.extend([
        {'parent': 'MOVERS',
        'tag': 'AddOrRemoveMatchCsts',
        'name': 'cstadd', 'cst_instruction': 'add_new'},
        {'parent': 'PROTOCOLS',
            'tag': 'Add',
            'mover_name': 'cstadd'}
        ])

    protocol_names:List[str] = list()
    for rd in reactant_dicts:
        cc: str = rd['chain']
        radius: float = rd['length']
        predock:str=f'predock_{cc.lower()}'
        dock:str=f'dock_{cc.lower()}'
        elements.extend([
        {
            'parent': 'MOVERS',
            'tag': 'Transform',
            'name': dock,
            'chain': f'{cc.upper()}',
            'box_size': f'{0.5 * radius:.3f}',
            'move_distance': '0.1',
            'angle': '360',
            'cycles': '1000',
            'repeats': '3',
            'temperature': '5',
            'grid_set': rd['grid_name']
        }])
        protocol_names.extend([dock])

    elements.extend(
        [{'parent': 'PROTOCOLS', 'tag': 'Add', 'mover_name': pn } for pn in protocol_names] +
        [{'parent': 'MOVERS',
            'tag': 'EnzRepackMinimize',
            'name':'cstopt',
            'minimize_rb':'0', 
            'cst_opt':'1',
            'cycles':'3',
            'scorefxn_repack': 'ligand_soft_rep',
            'scorefxn_minimize': 'ligand_soft_rep',
            'minimize_lig':'1',
            'minimize_sc':'1',
            'minimize_bb':'0',
            'cycles':'3',
            'min_in_stages':'0',},
        {'parent': 'PROTOCOLS',
            'tag': 'Add',
            'mover_name': 'cstopt'},
     ])
                

    interface.rosetta.write_script(fname, elements)

    _LOGGER.info(f"Saved new RosettaScripts .xml file at {fpath.absolute()}!")
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
        #"-keep_input_protonation_state", TODO(CJ): get rid of this for now; but want it back someday
        "-auto_setup_metals",
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
        "-enzdes:minimize_all_ligand_torsions 5.0",
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


def _dock_system(structure: Structure,
            constraints: List[RosettaCst],
            n_struct:int,
            param_files:List[str],
            charge_mapper:Dict[str,int],
            work_dir:str,
            rng_seed:int,
            use_cache:bool) -> pd.DataFrame:
    """

    Args:
        structure:
        constraints:
        n_struct:
        work_dir:
        rng_seed:
        use_cache:

    Returns:

    """
    
    if use_cache:
        csv_file = Path( f"{work_dir}/scores.csv")
        _LOGGER.info(f"Cache mode enabled. Looking for {csv_file} ...")
        if csv_file.exists():
            df = pd.read_csv(csv_file)
            df['selected'] = True
            _LOGGER.info(f"Found file {csv_file}. Checking for existence of {len(df)} output .pdbs")

            for i, row in df.iterrows():
                if not Path(row.description).exists():
                    _LOGGER.info(f"\t{row.description} is missing. Cannot use cached results!")
                    break
            else:
                if len(df) == n_struct:
                    _LOGGER.info("All cached .pdbs exist! Using cached RosettaLigand structures!")
                    return df
        else:
            _LOGGER.info(f"{csv_file} not found. Continuing with standard run")

    (start_pdb, cst_file) = _integrate_csts(structure, constraints, work_dir, use_cache)

    xml_file: str = _create_xml(structure, work_dir, use_cache)  #TODO(CJ): going to overhaul this

    options_file: str = _make_options_file(start_pdb, xml_file, param_files, work_dir, rng_seed, n_struct, use_cache, cst_file)

    opt_path = Path(options_file)

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

    return df, charge_mapper

def _evaluate_SASA(df: pd.DataFrame, max_sasa_ratio:float) -> None:
    """TODO(CJ)"""
    #TODO(CJ): put an error if the df is empty
    _LOGGER.info(f"Beginning ligand SASA evaluation. {df.selected.sum()} geometries still selected...")
    parser = PDBParser()
    stru:Structure = parser.get_structure(df.iloc[0].description)
   
    sele_str=str()
    ligand_keys = list()
    for residue in stru.residues:
        if not residue.is_ligand():
            continue
        ligand_keys.append((residue.chain.name, str(residue.idx)))
        sele_str += f"( chain {ligand_keys[-1][0]} and resi {ligand_keys[-1][1]}) or "

    session = interface.pymol.new_session()
    args = [('flag', 'ignore', 'none'), ('flag', 'ignore', 'solvent')]
    interface.pymol.general_cmd(session, args)
    good_sasa:List[bool]=list()
    for i,row in df.iterrows():
        if not row.selected:
            good_sasa.append( False )
        else:
            good = True
            args = [('delete', 'all'), ('load', row.description), ("get_sasa_relative", sele_str[:-3])]
            sasa_values = interface.pymol.general_cmd(session, args)[-1]
            for (_,_,sv_chain,sv_index),sasa_rel in sasa_values.items():
                if (sv_chain,sv_index) in ligand_keys:
                    if sasa_rel >= max_sasa_ratio:
                        good = False
            good_sasa.append(good)


    df['good_sasa'] = good_sasa
    df['selected'] = (df.selected) & (df.good_sasa)
    _LOGGER.info(f"Finished ligand SASA evaluation. {df.selected.sum()} geomtries have ligands with relative SASA <= {max_sasa_ratio:.3f}")


def _evaluate_csts(df: pd.DataFrame, csts: List[RosettaCst], cst_cutoff: int) -> None:
    """Evaluates geometries of the .pdb files in the supplied DataFrame versus the specified constraints.
    Constraints are evaluted in terms of how many tolerance units each specified angle, distance, etc. is
    different than the idealized value. The 'selected' column in the DataFrame will be updated by the 
    specified cst_cutoff value.

    Args:
        df: The DataFrame object taken from the 
        csts: A list() of RosettaCsts to evaluate against.
        cst_cutoff: The maximum amount of tolerance units allowed as an int().

    Returns:
        Nothing.
    """
    _LOGGER.info(f"Beginning RosettaCst evaluation. {df.selected.sum()} geometries still selected...")

    cst_diff = []
    for i, row in df.iterrows():
        total = 0.0
        for cst in csts:
            for tidx, tl in enumerate(cst.evaluate(row.description)):
                penalty = cst.constraints[tidx][3]
                if tl > 1.0:
                    total += (tl-1)*penalty
        cst_diff.append(total)

    df['cst_diff'] = cst_diff
    df['selected'] &= (df['cst_diff'] <= cst_cutoff)

    #exit( 0 )
    _LOGGER.info(
        f"Finished RosettaCst evaluation. {df.selected.sum()} geometries have constraint tolerances <= {cst_cutoff:.3f} tolerance units")


def _parameterize_system(stru:Structure, work_dir:str) -> Tuple[List[str], Dict[str, int]]:
    """Given the input Structure(), parameterize 
    Args:
        reactants:
        reactant_conformers:

    Returns:
        A Tuple() with the format (list(), dict()), holding the (.params file names, charge mapper ). 
    """
    _LOGGER.info("Beginning preparation of each reactant...")
    #TODO(CJ): error if not enough conformers... then make the conformers
    param_files: List[str] = list()
    charge_mapper: Dict[str, int] = dict()
   
    parser = Mol2Parser()
    for res in stru.residues:
        conformers:str=None
        if type(res) != Ligand:
            continue
        _LOGGER.info(f"Detected residue {res.name} in chain {res.parent.name}...")
        if res.n_conformers() == 1:
            _LOGGER.info(f"\tDetected no conformers!")
        else:
            _parser = PDBParser()
            
            orig:List[Tuple[float,float,float]] = [aa.coord for aa in res.atoms]

            content:List[str]=list()
            for conf_coords in res.conformer_coords:                
                for aidx, aa in enumerate(res.atoms):
                    aa.coord = conf_coords[aidx]
                content.extend(_parser.get_file_str(res, if_renumber=False).splitlines())

            for aidx, aa in enumerate(res.atoms):
                aa.coord = orig[aidx]
            
            conformers = f"{work_dir}/{res.name}_conformers.pdb"
            fs.write_lines(conformers, filter(lambda ll: not ll.startswith('TER'), content ) )

        temp_rct:str = f"{work_dir}/__temp_rct.mol2"
        parser.save_ligand(temp_rct, res)
        (pfile, charge) = _parameterize_reactant(
            temp_rct,
            res.name,
            conformers
        )
        _LOGGER.info(f"Information for reactant {res.name}:")
        _LOGGER.info(f"\tparam file: {pfile}")
        _LOGGER.info(f"\tcharge: {charge}")
        param_files.append(pfile)
        charge_mapper[res.name] = charge
        fs.safe_rm( temp_rct )

    _LOGGER.info("Finished reactant preparation!")

    return (param_files, charge_mapper)




def _ligand_chain_names(stru:Structure) -> List[str]:
    """TODO(CJ): 
    Args:

    Returns:
    """
    #TODO(CJ): change to reactant names
    _LOGGER.info(f"Quick analysis of reactant-enzyme system...")
    result: List[str] = list()

    _LOGGER.info("Detecting reactants:")
    for chain in stru.chains:
        for residue in chain.residues:
            if residue.is_ligand():
                _LOGGER.info(f"\treactant {residue.name} in chain {chain.name}")
                result.append(chain.name)

    return result

def _place_reactants(structure: Structure, reactants: List[str], pdb_code) -> bool:
    """TODO(CJ)"""

    _LOGGER.info("Beginning placement of reactants into apo-enzyme...")

    n_reactants: int = len(reactants)
    reactants_left: Set = set([Path(rr).stem for rr in reactants])

    for rr in structure.residues:
        if rr.name in reactants_left:
            reactants_left.remove(rr.name)

    if not reactants_left:
        _LOGGER.info("All reactants are present in complex! Continuing...")
        return

    if n_reactants == len(reactants_left):
        _LOGGER.info("Reactants not present in complex! Beginning placement strategy...")
        assert False, "Not implemented yet!!!!"
    else:
        _LOGGER.error("Some reactants are present in complex some are not. Must be all or nothing. Exiting...")
        exit(1)


def _evaluate_clashes(df: pd.DataFrame, clash_distance: float, clash_cutoff: int) -> None:
    """TODO(CJ)"""

    _LOGGER.info(f"Beginning clash evaluation. {df.selected.sum()} geometries still selected...")
    session = interface.pymol.new_session()
    clash_ct: List[int] = list()
    for i, row in df.iterrows():

        if not row.selected:
            clash_ct.append(-1)
            continue

        atoms: pd.DataFrame = interface.pymol.collect(session, row.description, "elem resn vdw x y z".split())
        atoms: pd.DataFrame = atoms[atoms.elem != 'H'].reset_index(drop=True)

        ligand_names = list(filter(lambda rn: rn[0] == 'L' and rn[1].isnumeric() and rn[2].isnumeric(), atoms.resn.unique()))

        other = atoms[~(atoms.resn.isin(ligand_names))].reset_index(drop=True)
        o_xyz = np.transpose(np.array([other.x.to_numpy(), other.y.to_numpy(), other.z.to_numpy()]))

        count: int = 0
        for resn in ligand_names:
            ligand = atoms[atoms.resn == resn].reset_index(drop=True)
            l_xyz = np.transpose(np.array([ligand.x.to_numpy(), ligand.y.to_numpy(), ligand.z.to_numpy()]))

            for ll in l_xyz:
                if np.min(np.sqrt(np.sum(np.power(o_xyz - ll, 2), axis=1))) <= clash_distance:
                    count += 1
        clash_ct.append(count)

    df['clash_ct'] = clash_ct
    df['selected'] = (df.selected) & (df.clash_ct <= clash_cutoff)

    _LOGGER.info(
        f"Finished clash evaluation. {df.selected.sum()} geometries have <= {clash_cutoff} of distance {clash_distance:.3f} angstroms")


def _create_binding_clusters(geometries: List[str], chain_names):
    _LOGGER.info("Creating binding clusters...")

    session = interface.pymol.new_session()

    rms_mapper = dict()
    sele_names = [Path(gg).stem for gg in geometries]

    args = []
    for gg in geometries:
        args.append(('load', gg))

    args.append(('remove', 'not ( ' + ' or '.join(map(lambda ll: f"chain {ll}", chain_names)) + ' ) '))
    args.append(('save', 'temp.pse'))

    interface.pymol.general_cmd(session, args)

    n_geo: int = len(geometries)
    chain_sele = "(" + " or ".join(map(lambda ll: f"chain {ll}", chain_names)) + ")"
    for gidx1 in range(n_geo):
        for gidx2 in range(gidx1 + 1, n_geo):
            g1, g2 = geometries[gidx1], geometries[gidx2]
            args = [('rms', f"{sele_names[gidx1]} and {chain_sele}", f"{sele_names[gidx2]} and {chain_sele}")]
            rms_value: float = interface.pymol.general_cmd(session, args)[-1]
            rms_mapper[(g1, g2)] = rms_value
            rms_mapper[(g2, g1)] = rms_value
    interface.pymol.general_cmd(session, [('delete', 'all')])

    rms_clusters = list()

    for gg in geometries:
        placed = False
        for grp in rms_clusters:
            for member in grp:
                if rms_mapper[(gg, member)] <= 2.0:  #TODO(CJ): paramterize this
                    grp.add(gg)
                    placed = True
                if placed:
                    break
            if placed:
                break
        if not placed:
            rms_clusters.append({gg})

    _LOGGER.info(f"Created {len(rms_clusters)} clusters!")

    return rms_clusters


def _system_charge(df: pd.DataFrame, charge_mapper=None) -> int:
    """TODO(CJ) """
    _LOGGER.info("Beginning enzyme system charge deduction...")
    charge = 0

    residues = sorted(list(set(zip(df.chain, df.resi, df.resn))))
    _LOGGER.info(f"Found {len(residues)} residues in enzyme system!")

    for (chain, resi, resn) in residues:
        atoms = df[(df.chain == chain) & (df.resi == resi)].name.to_list()

        charge_str = str()
        residue_charge: int = 0
        if resn in chem.AA_CHARGE_MAPPER:
            residue_charge = chem.AA_CHARGE_MAPPER[resn]
            charge_str = f"Detected {resn} with charge: {residue_charge}"
        elif resn in charge_mapper:
            residue_charge = charge_mapper[resn]
            charge_str = f"Detected {resn} with charge: {residue_charge}"
        elif resn.upper() == "MG":  #TODO(CJ): add more charges
            residue_charge = 2
            charge_str = f"Detected metal ion {resn} with charge: {residue_charge}"
        else:
            _LOGGER.error(f"Unable to assign charge to {resn}. Exiting...")
            exit(1)

        # correcting titratable residues
        special_titration: bool = False
        if resn == "GLU":
            if "HE2" in atoms:
                residue_charge += 1
                special_titration = True
            if "HE1" in atoms:
                residue_charge += 1
                special_titration = True

        if resn == "LYS":
            if "3HZ" not in atoms:
                residue_charge -= 1
                special_titration = True

        if resn == "ASP":
            if "HD2" in atoms:
                residue_charge += 1
                special_titration = True
            if "HD1" in atoms:
                residue_charge += 1
                special_titration = True

        if resn == "TYR":
            if "HH" not in atoms:
                residue_charge -= 1
                special_titration = True

        if resn == "HIS":
            if "HD1" in atoms and "HE2" in atoms:
                residue_charge += 1
                special_titration = True

        if special_titration:
            charge_str = f"Detected special protonation state of {resn}. New charge: {residue_charge}"

        charge += residue_charge
        _LOGGER.info(charge_str)

    _LOGGER.info(f"Total charge of enzyme system: {charge}")

    return charge


def _define_binding_pockets(start_pdb: str, distance_cutoff, charge_mapper=None):
    """TODO"""
    #TODO(CJ): make this work for multiple reactants
    session = interface.pymol.new_session()
    res_names: List[str] = interface.pymol.collect(session, start_pdb, "resn".split()).resn.unique()
    sele_names: Dict[str, str] = {}
    binding_ddg: Dict[str, List[float]] = {}
    charges: Dict[str, float] = {}

    binding_pockets = list()

    reactants = list()
    for rn in res_names:
        if rn.upper() in chem.METAL_CENTER_MAP or rn.upper() in chem.THREE_LETTER_AA_MAPPER:
            continue
        reactants.append(rn)

    _LOGGER.info(f"Found {len(rn)} reactants: {', '.join(reactants)}")
    _LOGGER.info("Analyzing binding pockets of each reactant...")

    reactant_sele = " or ".join(map(lambda ll: f"resn {ll}", reactants))
    for rn in reactants:

        binding_pocket = dict()
        _LOGGER.info(f"Analyzing reactant {rn}...")
        atoms: pd.DataFrame = interface.pymol.collect(
            session,
            start_pdb,
            "resi chain resn name".split(),
            sele=f"(byres all within {distance_cutoff:.2f} of resn {rn}) and not ({reactant_sele})")

        binding_pocket["receptor_charge"] = _system_charge(atoms, charge_mapper)
        binding_pocket["probe_charge"] = charge_mapper[rn]
        binding_pocket["sele"] = " or ".join(map(lambda pr: f"(resi {pr[0]} and chain {pr[1]})", set(zip(atoms.resi, atoms.chain))))
        binding_pocket["resn"] = rn

        _LOGGER.info(f"\tprobe charge: {binding_pocket['probe_charge']}")
        _LOGGER.info(f"\treceptor charge {binding_pocket['receptor_charge']}")
        _LOGGER.info(f"\tnumber of residues in receptor: {binding_pocket['sele'].count(' or ')+1}")

        binding_pockets.append(binding_pocket)

    return binding_pockets


def _evaluate_binding(df: pd.DataFrame,
                      start_pdb: str,
                      binding_cutoff: int,
                      chain_names,
                      distance_cutoff: float = 4.0,
                      use_rms: bool = True,
                      charge_mapper=None) -> None:
    """ """

    _LOGGER.info(f"Beginning binding evaluation. {df.selected.sum()} geometries still selected...")

    rms_clusters = _create_binding_clusters(df[df.selected].description.to_list(), chain_names)

    binding_pockets: List[Set[str]] = _define_binding_pockets(start_pdb, distance_cutoff, charge_mapper)

    clusters = dict()
    cluster_mapper = dict()
    resnames = set()
    for cidx, rcluster in enumerate(rms_clusters):
        cluster_name: str = f"cluster_{cidx:03d}"
        cluster = {"name": cluster_name}
        _LOGGER.info(f"Analyzing cluster {cidx+1} of {len(rms_clusters)}..")
        molfile: str = list(rcluster)[0]
        for rc in rcluster:
            cluster_mapper[rc] = cluster_name

        _LOGGER.info(f"Using file {molfile}")

        for bp in binding_pockets:
            #TODO(CJ): check if xtb failed for SCF iteration reasons and retry if so
            _LOGGER.info(f"Calculating binding energy for {bp['resn']}...")
            be: float = binding_energy(molfile,
                                       f"resn {bp['resn']}",
                                       bp['sele'],
                                       bp['probe_charge'],
                                       bp['receptor_charge'],
                                       work_dir=config["system.SCRATCH_DIR"])
            _LOGGER.info(f"Found binding energy of {be:.3f} hartrees")
            cluster[bp['resn']] = be
            resnames.add(bp['resn'])

        clusters[cluster_name] = cluster

    new_data = defaultdict(list)

    for i, row in df.iterrows():
        cluster_name: str = cluster_mapper.get(row.description, None)
        new_data['cluster'].append(cluster_mapper.get(row.description, None))
        for rn in resnames:
            temp = clusters.get(cluster_name, None)
            if temp is not None:
                temp = temp.get(rn, None)
            new_data[f"{rn}_ddg"].append(temp)

    binding_log_sum = np.zeros(len(df))
    for col_name, col_values in new_data.items():
        df[col_name] = col_values

        if col_name != 'cluster':
            binding_log_sum += np.exp(df[col_name].to_numpy())

    df['binding_ddg_log_sum'] = binding_log_sum

    bls = df.binding_ddg_log_sum.to_numpy()
    bls_is_nan = np.isnan(bls)
    bls_cutoff = np.percentile(bls[~bls_is_nan], binding_cutoff)
    df['selected'] = (df.selected) & (~bls_is_nan) & (bls <= bls_cutoff)

    _LOGGER.info(f"Finished binding evaluation. {df.selected.sum()} geometries have been selected.")


def _define_active_site(structure: Structure, distance_cutoff, charge_mapper=None):
    """TODO(CJ)"""

    parser = PDBParser()
    start_pdb:str=f"{config['system.SCRATCH_DIR']}/__temp_as.pdb"
    parser.save_structure(start_pdb, structure)

    session = interface.pymol.new_session()
    res_names: List[str] = interface.pymol.collect(session, start_pdb, "resn".split()).resn.unique()
    sele_names: Dict[str, str] = {}
    binding_ddg: Dict[str, List[float]] = {}
    charges: Dict[str, float] = {}

    binding_pockets = list()

    _LOGGER.info("Analyzing enzyme active site...")
    reactants = list()
    sele = []
    for rn in res_names:
        if rn.upper() in chem.METAL_CENTER_MAP or rn.upper() in chem.THREE_LETTER_AA_MAPPER:
            continue
        reactants.append(rn)
        sele.append(f"((byres all within {distance_cutoff:.2f} of resn {rn}) or (metals within  {2*distance_cutoff:.2f} of resn {rn}))")
    sele_str = " or ".join(sele)
    atoms: pd.DataFrame = interface.pymol.collect(session, start_pdb, "resi chain resn name".split(), sele=sele_str)
    residues = sorted(list(set(zip(atoms.chain, atoms.resi))))
    charge = _system_charge(atoms, charge_mapper)
    _LOGGER.info(f"Found {len(residues)} residues within {distance_cutoff} angstroms of reactants. Active site info:")
    _LOGGER.info(f"\tresidues: {len(residues)}")
    _LOGGER.info(f"\tcharge: {charge}")

    fs.safe_rm(start_pdb)

    sele_str = " or ".join(map(lambda pr: f"(chain {pr[0]} and resi {pr[1]})", residues))

    return {"sele": sele_str, "charge": charge}


def _qm_minimization(df:pd.DataFrame, charge_mapper, cluster_distance:float, constraints=None, work_dir:str=None) -> str:
    """TODO(CJ)
    
    Args:
        df:
        charge_mapper:
        cluster_distance:
        constraints:
        work_dir:

    Returns:

    """
    def __map(asite, enzyme):
        mapper = dict()

        atom_mapper:Dict = dict(zip(
            zip(enzyme.chain, enzyme.resn, enzyme.resi, enzyme['name']),enzyme['rank']
        ))

        for aidx, arow in asite.iterrows():
            a_key = (arow.chain, arow.resn, arow.resi, arow['name'])

            result = atom_mapper.get(a_key, None)
            mapper[(arow.chain, arow.resn, int(arow.resi), arow['name'])] = arow['rank']

            if result is None:
                continue

            #TODO(CJ): map more stuff
            mapper[arow['rank']] = result

        return mapper

    if work_dir is None:
        work_dir = config['system.SCRATCH_DIR']
    _LOGGER.info(df.sort_values(by='qm_energy'))
    infile:str=df[df.selected].sort_values(by='qm_energy').description.to_list()[0]
    
    session = interface.pymol.new_session()
    non_aa:pd.DataFrame = interface.pymol.collect(session, infile, 'chain resn resi'.split(), 'not polymer.protein')
    
    sele:str = " or ".join(
        map(lambda pr: f"( byres all within {cluster_distance} of  (chain {pr[0]} and resn {pr[1]} and resi {pr[2]}) )",
        set(zip(non_aa.chain, non_aa.resn, non_aa.resi))
    ))
    
    full_system:pd.DataFrame = interface.pymol.collect(session, infile, 'chain resn resi name rank'.split())
    
    temp_pdb = f"{work_dir}/__qm_min_temp.pdb"
    temp_mol = f"{work_dir}/__qm_min_temp.mol"
    
    temp_pdb = interface.pymol.create_cluster(session, infile, sele_str=sele, outfile=temp_pdb, cap_strategy='CH3')
    active_site:pd.DataFrame = interface.pymol.collect(session, temp_pdb, 'chain resn resi name rank'.split(), sele)

    atom_mapper:Dict = __map(active_site, full_system)
    charge:int = _system_charge(active_site, charge_mapper)

    temp_df = interface.pymol.collect(interface.pymol.new_session(), temp_pdb, 'chain resn resi name rank'.split())
    
    freeze_atoms = interface.pymol.collect(interface.pymol.new_session(), temp_pdb, 'rank'.split(), 'bb.')['rank'].to_numpy() + 1

    temp_sdf = interface.pymol.create_cluster(session, infile, sele_str=sele, outfile=temp_mol, cap_strategy='CH3')

    mapped_constraints:List = list()
    if constraints is not None:
        for cst in constraints:
            for raw_cst in cst.get_constraints():
                mapped_atoms = list()
                for ra in raw_cst['atoms']:
                    mapped_atom = atom_mapper.get(ra, None)
                    if mapped_atom is None:
                        continue
                    mapped_atoms.append(mapped_atom + 1)
                if len(mapped_atoms) != len(raw_cst['atoms']):
                    continue
                mapped_constraints.append((raw_cst['generic_type'], raw_cst['ideal_value'], mapped_atoms))


    (optimized_cluster, energy) = interface.xtb.geo_opt(
        temp_mol,
        charge=charge,
        freeze_atoms=freeze_atoms,
        constraints=mapped_constraints
        )

    optimized_df = interface.pymol.collect(interface.pymol.new_session(), optimized_cluster, 'x y z rank'.split())

    args = [('load', infile)]
    for i, row in optimized_df.iterrows():
        mapped_index = atom_mapper.get(row['rank'])
        if mapped_index is None:
            continue

        args.append((
            'alter_state', 1, f'(rank {mapped_index})',f'(x,y,z)=({row["x"]},{row["y"]},{row["z"]})'
        ))
        
    temp_path = Path(infile)
    outfile = str(temp_path.parent / f"{temp_path.stem}_optimized.pdb")
    args.append(('save', outfile))

    session = interface.pymol.new_session()
    interface.pymol.general_cmd(session, args)
    return outfile

def _evaluate_qm(df: pd.DataFrame, structure: Structure, charge_mapper, cluster_cutoff: float) -> None:
    """TODO(CJ) """
    # steps
    # 1. get system charge
    # 2. run through each cluster and do it
    _LOGGER.info(f"Beginning qm energy evaluation. {df.selected.sum()} geometries still selected...")
    as_info = _define_active_site(structure, cluster_cutoff, charge_mapper)
    session = interface.pymol.new_session()

    qm_energy = []

    for i, row in df.iterrows():
        if not row.selected:
            qm_energy.append(None)
            continue

        interface.pymol.general_cmd(session, [('delete', 'all')])
        interface.pymol.create_cluster(session, row.description, as_info['sele'], outfile='temp.xyz', cap_strategy='CH3')
        qm_energy.append(interface.xtb.single_point('temp.xyz', charge=as_info['charge']))

        fs.safe_rm('temp.xyz')
        _LOGGER.info(row.description)

    df['qm_energy'] = qm_energy
    _LOGGER.info("Finished qm energy evaluation!")


def _log_settings(
    constraints:List[RosettaCst],
    n_struct: int,
    cst_energy:float,
    clash_distance:float,
    clash_cutoff:int,
    max_sasa_ratio:float,
    cluster_distance:float,
    use_cache:bool,
    rng_seed:int,
    work_dir:str,
    save_work_dir:bool
):
    """Logs settings for the current reactive docking run. Done at the INFO level.

    Args:

    Returns:
        Nothing.
    """
    _LOGGER.info("Beginning EnzyRCD Reactive docking run! Below are the run settings and characteristics:")

    if constraints is not None:
        _LOGGER.info(f"\t{len(constraints)} RosettaConstraints")
    else:
        _LOGGER.info("\t0 RosettaConstraints")

    _LOGGER.info(f"\t{n_struct=}")
    _LOGGER.info(f"\t{cst_energy=:.3f} rosetta energy units")
    _LOGGER.info(f"\t{clash_distance=:.3f} angstroms")
    _LOGGER.info(f"\t{clash_cutoff=} clashes")
    _LOGGER.info(f"\t{max_sasa_ratio=:.3f} maximum allowed SASA ratio")
    _LOGGER.info(f"\t{cluster_distance=:.3f} angstroms")

    _LOGGER.info(f"\t{use_cache=}")
    _LOGGER.info(f"\t{rng_seed=}")
    _LOGGER.info(f"\t{work_dir=}")
    _LOGGER.info(f"\t{save_work_dir=}")
