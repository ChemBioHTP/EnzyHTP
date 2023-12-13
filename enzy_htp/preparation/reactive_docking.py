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
import enzy_htp.structure.structure_operation as stru_oper
from enzy_htp.structure.structure_constraint import StructureConstraint

import enzy_htp.chemical as chem
from enzy_htp import mutation as mm
from enzy_htp.structure import PDBParser, Mol2Parser, Structure, Ligand
from enzy_htp.core import file_system as fs

#TODO(CJ): I should go through and use existing values if cache mode is selected
#:TODO(CJ): need to be able to specify the size of the grid

def dock_reactants(structure: Structure,
                   constraints: List[StructureConstraint] = None,
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

    interface.rosetta.rename_atoms(structure)

    _log_settings(constraints, n_struct, cst_energy, clash_distance, clash_cutoff, max_sasa_ratio,
                cluster_distance, use_cache, rng_seed, work_dir, save_work_dir)
    
    (param_files, charge_mapper) = _parameterize_system(structure, work_dir)

    df = _dock_system(structure, constraints, n_struct, param_files, charge_mapper, work_dir, rng_seed, use_cache )
    
    _evaluate_SASA(df, max_sasa_ratio)#TODO(CJ): put this back in 

    _evaluate_clashes(df, clash_distance, clash_cutoff)

    _evaluate_csts(df, constraints, cst_energy)
    
    _evaluate_qm(df, structure, charge_mapper, cluster_distance)

    final_geometry:Structure = _qm_minimization(df, charge_mapper, cluster_distance, constraints)
   
    if not save_work_dir:
        _LOGGER.info(f"save_work_dir set to False! Deleting {work_dir}")
        fs.safe_rmdir( work_dir )

    return final_geometry 



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
            'grid_set': rd['grid_name'],
            'use_constraints':'true',
            'cst_fa_file':'./dock.cst'
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
        "-keep_input_protonation_state", #TODO(CJ): get rid of this for now; but want it back someday
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
            constraints: List[StructureConstraint],
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

    (start_pdb, cst_file) = interface.rosetta.integrate_constraints(structure, constraints, work_dir)

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

    return df

def _evaluate_SASA(df: pd.DataFrame, max_sasa_ratio:float) -> None:
    """Evaluates and filters geometries using SASA ratio of ligands. Updates the 'selected'
    column of the inputted DataFrame and sets the column to False if the SASA ratio of the ligand is above
    the specified cutoff. Exits if no rows are selected at this point.
    
    Args:
        df: pandas DataFrame with geometries derived from _dock_system().
        max_sasa_ratio: The SASA cutoff for the ligand.

    Returns:
        Nothing.
    """
    if not df.selected.sum():
        _LOGGER.error("No geometries still selected! Exiting...")
        exit( 1 )

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


def _evaluate_csts(df: pd.DataFrame, csts: List[StructureConstraint], cst_cutoff: int) -> None:
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
    _parser = PDBParser()
    cst_diff = []
    for i, row in df.iterrows():
        if not row.selected:
            cst_diff.append( None )
            continue

        total:float = 0.0
        stru:Structure = _parser.get_structure(row.description)
        for cst in csts:
            cst.change_topology(stru)
            total += cst.score_energy()
        cst_diff.append(total)

    df['cst_diff'] = cst_diff
    df['selected'] &= (df['cst_diff'] <= cst_cutoff)

    _LOGGER.info(
        f"Finished RosettaCst evaluation. {df.selected.sum()} geometries have constraint tolerances <= {cst_cutoff:.3f} tolerance units")
    

def _parameterize_system(stru:Structure, work_dir:str) -> Tuple[List[str], Dict[str, int]]:
    """Given the input Structure(), parameterize everything needed to perform RosettaLigand docking.
    
    Args:
        stru: The Structure() to parameterize.
        work_dir: Where temporary files will be saved.

    Returns:
        A Tuple() with the format (list(), dict()), holding the (.params file names, charge mapper ). 
    """
    _LOGGER.info("Beginning preparation of each reactant...")
    param_files: List[str] = list()
    charge_mapper: Dict[str, int] = dict()
   
    parser = Mol2Parser()
    for res in stru.residues:
        conformers:str=None
        if type(res) != Ligand:
            continue
        
        _LOGGER.info(f"Detected residue {res.name} in chain {res.parent.name}...")

        ligand_charge:int = charge_mapper.get(res.name, None)

        if ligand_charge is None:
            ligand_charge = interface.bcl.calculate_formal_charge( res )
            charge_mapper[res.name] = ligand_charge
        _LOGGER.info(charge_mapper)

        param_file:str = interface.rosetta.parameterize_ligand(res, charge=ligand_charge, work_dir=work_dir)
        
        _LOGGER.info(f"Information for reactant {res.name}:")
        _LOGGER.info(f"\tparam file: {param_file}")
        _LOGGER.info(f"\tcharge: {ligand_charge}")
        _LOGGER.info(f"\tnum conformers: {res.n_conformers()}")
        
        param_files.append(param_file)
        charge_mapper[res.name] = ligand_charge
    _LOGGER.info("Finished reactant preparation!")

    return (param_files, charge_mapper)


def _evaluate_clashes(df: pd.DataFrame, clash_distance: float, clash_cutoff: int) -> None:
    """Evaluates and filters geometries using clash count between each ligand and other ligands.
    Updates the 'selected' column of the inputted DataFrame and sets the column to False if the 
    ligand has more than the allowed number of clashes between heavy atoms. Exits if no rows are selected at this point.
    
    Args:
        df: pandas DataFrame with geometries derived from _dock_system().
        clash_distance: Allowed clash distance radius in Angstroms.
        clash_cutoff: Allowed number of clashes.

    Returns:
        Nothing.         
    """
    #TODO(CJ): there are problems with this since we are no longer using LXX format ligands 
    if not df.selected.sum():
        _LOGGER.error("No geometries still selected! Exiting...")
        exit( 1 )

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
    if work_dir is None:
        work_dir = config['system.SCRATCH_DIR']

    if not df.selected.sum():
        _LOGGER.error("No geometries still selected! Exiting...")
        exit( 1 )
    
    #TODO(CJ): add the xtb.geo_opt() here

    infile:str=df[df.selected].sort_values(by='qm_energy').description.to_list()[0]
   
    _parser = PDBParser()
    final_stru = _parser.get_structure(infile)
    as_info = _define_active_site(final_stru, cluster_distance, charge_mapper)
    
    interface.xtb.geo_opt(final_stru, charge=as_info['charge'], constraints=constraints, sele_str=as_info['sele'], work_dir=work_dir)

    return final_stru

def _evaluate_qm(df: pd.DataFrame, structure: Structure, charge_mapper, cluster_cutoff: float) -> None:
    """TODO(CJ)

    Args:
        df:
        structure:
        charge_mapper:
        cluster_cutoff:

    Returns:
        Nothing.        

    """
    # steps
    # 1. get system charge
    # 2. run through each cluster and do it
    _LOGGER.info(f"Beginning qm energy evaluation. {df.selected.sum()} geometries still selected...")
    as_info:Dict = _define_active_site(structure, cluster_cutoff, charge_mapper)

    qm_energy = []

    _parser = PDBParser()

    for i, row in df.iterrows():
        if not row.selected:
            qm_energy.append(None)
            continue

        energy:float = None
        _df_stru = _parser.get_structure( row.description )
        try:
            energy = interface.xtb.single_point(_df_stru,charge=as_info['charge'],sele_str=as_info['sele'])
        except:     
            _LOGGER.warn("Failed single point energy calculation attempt! Continuing...")
            #TODO(CJ): add some file cleanup

        qm_energy.append(energy)

    df['qm_energy'] = qm_energy
    _LOGGER.info("Finished qm energy evaluation!")


def _log_settings(
    constraints:List[StructureConstraint],
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
        constraints: A List() of RosettaCst's defining the system at hand.
        n_struct: How many geometries to sample as an int(). 
        cst_energy: Rosetta Energy Unit cutoff for filtering geometries.
        clash_distance: Heavy atom radius for clash counting as a float().
        clash_cutoff: The number of clashes before a geometry is filtered as an int().
        max_sasa_ratio: Maximum allowed SASA ratio for ligands in a geometry as a float().
        use_cache: Should existing results be used when possible?
        rng_seed: The seed for random number generation as an int(). 
        work_dir: Directory where temporary files are written. 
        save_work_dir: Should temporary files be deleted after reactive docking is run.?

    Returns:
        Nothing.
    """
    _LOGGER.info("Beginning EnzyRCD Reactive docking run! Below are the run settings and characteristics:")

    if constraints is not None:
        _LOGGER.info(f"\t{len(constraints)} StructureConstraints")
    else:
        _LOGGER.info("\t0 StructureConstraints")

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
