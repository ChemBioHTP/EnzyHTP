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

from enzy_htp import interface, config, _LOGGER
import enzy_htp.structure.structure_operation as stru_oper
from enzy_htp.structure.structure_constraint import StructureConstraint

import enzy_htp.chemical as chem
from enzy_htp import mutation as mm
from enzy_htp.structure import PDBParser, Mol2Parser, Structure, Ligand
from enzy_htp.core import file_system as fs
from enzy_htp.quantum import single_point

#TODO(CJ): I should go through and use existing values if cache mode is selected
#:TODO(CJ): need to be able to specify the size of the grid

def dock_reactants(structure: Structure,
                   constraints: List[StructureConstraint] = None,
                   n_struct: int = 100,
                   cst_energy: float = None,
                   use_qm: bool = True,
                   freeze_alphafill:bool = True,
                   clash_distance: float = 2.0,
                   clash_cutoff: int = 3,
                   max_sasa_ratio:float = 0.60,
                   cluster_distance: float = 2.0,
                   rng_seed: int = 1996,
                   work_dir: str = None,
                   save_work_dir: bool = True,
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
        use_qm:
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
    
    (param_files, charge_mapper) = _parameterize_system(structure, work_dir)
    
    geometry_df = pd.read_csv("scratch/scores.csv")

    if False:
        geometry_df = generate_geometries(structure, constraints, n_struct, clash_cutoff, param_files, charge_mapper, freeze_alphafill, rng_seed, work_dir )

    select_geometry(structure, geometry_df, max_sasa_ratio, constraints, cst_energy, charge_mapper, cluster_distance)

    if use_qm:
        qm_minimization(structure, charge_mapper, cluster_distance, constraints)

    mm_minimization(structure, cluster_distance, constraints, param_files, rng_seed, work_dir)

    if use_qm:
        qm_minimization(structure, charge_mapper, cluster_distance, [])

    if not save_work_dir:
        _LOGGER.info(f"save_work_dir set to False! Deleting {work_dir}")
        fs.safe_rmdir( work_dir )

    return structure 


def _create_xml(stru:Structure, cst_file:str, freeze_alphafill:bool, clash_cutoff:int, work_dir: str) -> str:
    """

    Args:
        

    Returns:
        The relative filepath of the RosettaScripts .xml file.
    """
    fname: str = f"{work_dir}/__script.xml"
    
    fpath = Path(fname)
#yapf: disable
    elements: List[Dict[str,str]] = [
        {'parent': 'SCOREFXNS', 'tag': 'ScoreFunction', 'name': 'ligand_soft_rep', 'weights': 'ligand_soft_rep'},
        {'parent': 'SCOREFXNS', 'tag': 'ScoreFunction', 'name': 'hard_rep', 'weights': 'ligand'},
        {'parent': 'SCOREFXNS.ScoreFunction', 'tag': 'Reweight', 'scoretype': 'coordinate_constraint',  'weight': '1.0'},
        {'parent': 'SCOREFXNS.ScoreFunction', 'tag': 'Reweight', 'scoretype': 'atom_pair_constraint',   'weight': '1.0'},
        {'parent': 'SCOREFXNS.ScoreFunction', 'tag': 'Reweight', 'scoretype': 'angle_constraint',       'weight': '1.0'},
        {'parent': 'SCOREFXNS.ScoreFunction', 'tag': 'Reweight', 'scoretype': 'dihedral_constraint',    'weight': '1.0'},
        {'parent': 'SCOREFXNS.ScoreFunction', 'tag': 'Reweight', 'scoretype': 'chainbreak',             'weight': '1.0'},
    ]
#yapf: enable
    
    as_names:List[str] = list()
    for res in stru.residues:
        if not res.is_ligand():
            continue

        chain_name:str = res.parent.name        
        rname:str=f"rs_{chain_name.lower()}"
        transform_name:str=f"dock_{chain_name.lower()}"
        elements.extend([
            {'parent': 'RESIDUE_SELECTORS', 'tag':'Index', 'name':rname, 'resnums':f"{res.idx}{res.parent.name}"},
            {'parent': 'RESIDUE_SELECTORS', 'tag':'CloseContact', 'name':f"{rname}_as", 'residue_selector':rname},
        ])

        as_names.append( f"{rname}_as" )

    elements.extend([
        { 'parent': 'RESIDUE_SELECTORS', 'tag': 'Or', 'name':'as_selector', 'selectors':','.join(as_names)},
        { 'parent': 'RESIDUE_SELECTORS', 'tag': 'Not', 'name':'not_as_selector', 'selector':'as_selector'},
        {'parent':'MOVERS', 'tag':'FastRelax', 'name':'frelax', 'scorefxn':'hard_rep', 'cst_file':f"{Path(cst_file).absolute()}", 'child_nodes': [
                {'parent':'FastRelax', 'tag':'MoveMap', 'name':'full_enzyme', 'bb':'true', 'chi':'true', 'jump':'false','child_nodes':[
                        {'parent':'MoveMap','tag':'ResidueSelector', 'selector':'as_selector', 'bb':'true', 'chi':'true', 'bondangle':'true'},
                        {'parent':'MoveMap', 'tag':'ResidueSelector', 'selector':'not_as_selector', 'bb':'false', 'chi':'true', 'bondangle':'false'},]}]
        },
        {'parent':'PROTOCOLS', 'tag':'Add', 'mover_name':'frelax'}
    ])

    interface.rosetta.write_script(fname, elements)
    _LOGGER.info(f"Saved new RosettaScripts .xml file at {fpath.absolute()}!")
    return fname


def _make_options_file(pdb_file: str,
                      xml_file: str,
                      param_files: List[str],
                      work_dir: str,
                      rng_seed: int,
                      n_struct: int,) -> str:
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
    return fname #TODO(CJ): fix this!
    score_file: str = f"{work_dir}/complexes/score.sc"

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


def generate_geometries(structure: Structure,
            constraints: List[StructureConstraint],
            n_struct:int,
            clash_cutoff:int,
            param_files:List[str],
            charge_mapper:Dict[str,int],
            freeze_alphafill:bool,
            rng_seed:int,
            work_dir:str,
            ) -> pd.DataFrame:
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
    start_pdb:str = f"{work_dir}/start.pdb"
    parser = PDBParser()
    parser.save_structure(start_pdb, structure)

    cst_file:str = interface.rosetta.write_constraint_file(structure, constraints, work_dir) #TODO(CJ): look at this; wrong constraint types!!

    xml_file:str = _create_xml(structure, cst_file, freeze_alphafill, clash_cutoff, work_dir)  #TODO(CJ): going to overhaul this

    options_file: str = _make_options_file(start_pdb,
                        xml_file, param_files, work_dir, rng_seed, n_struct)

    opt_path = Path(options_file)

    _LOGGER.info("Beginning RosettaLigand geometry sampling step...")
    
    start_dir: str = os.getcwd()
    
    os.chdir(str(opt_path.parent))
    
    try:
        interface.rosetta.run_rosetta_scripts([f"@{opt_path.name}"])
    except:
        pass

    os.chdir(start_dir)

    assert Path(str((opt_path.parent / "complexes/score.sc").absolute())).exists()

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

        param_file:str = interface.rosetta.parameterize_ligand(res, charge=ligand_charge, work_dir=work_dir)
        
        _LOGGER.info(f"Information for reactant {res.name}:")
        _LOGGER.info(f"\tparam file: {param_file}")
        _LOGGER.info(f"\tcharge: {ligand_charge}")
        _LOGGER.info(f"\tnum conformers: {res.n_conformers()}")
        
        param_files.append(param_file)
        charge_mapper[res.name] = ligand_charge
    _LOGGER.info("Finished reactant preparation!")

    return (param_files, charge_mapper)


def _system_charge(df: pd.DataFrame, charge_mapper:Dict[str,int]=None) -> int:
    """Given a pd.DataFrame with the columns chain, resi, resn, and name, find the charge. Needs a charge
    mapper for non-standard residues. If the charge mapper is missing, may be inaccurate.

    Args:
        df: The dataframe containing the chain, resi, resn, and name attributes for each atom in the selection.
        charge_mapper: Optional. A dict() that maps three letter residues to their charges.

    Returns:
        The system's charge as an int().
    """
    _LOGGER.info("Beginning enzyme system charge deduction...")
    charge = 0

    residues = sorted(list(set(zip(df.chain, df.resi, df.resn))))
    _LOGGER.info(f"Found {len(residues)} residues in enzyme system!")

    for (chain, resi, resn) in residues:
        atoms = df[(df.chain == chain) & (df.resi == resi)]['name'].to_list()

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

def _define_active_site(structure: Structure, distance_cutoff:float, charge_mapper:Dict[str,int]=None) -> Dict:
    """
    
    Args:
        
    Returns:

    """

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

def qm_minimization(structure:Structure, charge_mapper, cluster_distance:float, constraints:List[StructureConstraint]=None, work_dir:str=None):
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

    #TODO(CJ): check if the constraints belong to the Structure 
    as_info = _define_active_site(structure, cluster_distance, charge_mapper)
    
    interface.xtb.geo_opt(structure, charge=as_info['charge'], constraints=constraints, sele_str=as_info['sele'], work_dir=work_dir)

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
        for res in structure.residues:
            if res.is_canonical():
                continue
            _df_stru.get(res.key_str).net_charge = res.net_charge
            _df_stru.get(res.key_str).multiplicity = res.multiplicity
            for atom in _df_stru.get(res.key_str):
                atom.charge = 0.0
        lot = chem.QMLevelOfTheory(basis_set='',method='GFN2', solvent='', solv_method='')
        es = single_point(
            _df_stru,
            engine='xtb',
            region_methods=[lot],
            parallel_method=None,
            regions=[as_info['sele']])
        qm_energy.append(es[0].energy_0)

    df['qm_energy'] = qm_energy
    _LOGGER.info("Finished qm energy evaluation!")

def select_geometry( structure, df, max_sasa_ratio, constraints, cst_energy, charge_mapper, cluster_distance) -> None:
    """TODO(CJ)"""
    #_evaluate_SASA(df, max_sasa_ratio)#TODO(CJ): put this back in 

    #_evaluate_csts(df, constraints, cst_energy)
    #TODO(CJ): add use_qm flag 
    _evaluate_qm(df, structure, charge_mapper, cluster_distance)

    infile:str=df[df.selected].sort_values(by='qm_energy').description.to_list()[0]

    _parser = PDBParser()
    ref_stru = _parser.get_structure(infile)
    stru_oper.update_residues(structure, ref_stru)

    for cst in constraints:
        cst.change_topology(structure)

def mm_minimization(structure:Structure, cluster_distance:float, constraints:List[StructureConstraint], param_files:List[str], rng_seed:int, work_dir:str) -> None:
    """TODO(CJ)
    """

    cst_file:str = interface.rosetta.write_constraint_file(structure, constraints, work_dir) #TODO(CJ): look at this; wrong constraint types!!
    xml_file:str = f"{work_dir}/__script.xml"
    pdb_file:str = f"{work_dir}/__temp_mm_min.pdb"
    option_file:str = f"{work_dir}/__mm_min_options.txt"
#yapf: disable
    elements: List[Dict[str,str]] = [
        {'parent': 'SCOREFXNS', 'tag': 'ScoreFunction', 'name': 'ligand_soft_rep', 'weights': 'ligand_soft_rep'},
        {'parent': 'SCOREFXNS', 'tag': 'ScoreFunction', 'name': 'hard_rep', 'weights': 'ligand'},
        {'parent': 'SCOREFXNS.ScoreFunction', 'tag': 'Reweight', 'scoretype': 'coordinate_constraint',  'weight': '1.0'},
        {'parent': 'SCOREFXNS.ScoreFunction', 'tag': 'Reweight', 'scoretype': 'atom_pair_constraint',   'weight': '1.0'},
        {'parent': 'SCOREFXNS.ScoreFunction', 'tag': 'Reweight', 'scoretype': 'angle_constraint',       'weight': '1.0'},
        {'parent': 'SCOREFXNS.ScoreFunction', 'tag': 'Reweight', 'scoretype': 'dihedral_constraint',    'weight': '1.0'},
        {'parent': 'SCOREFXNS.ScoreFunction', 'tag': 'Reweight', 'scoretype': 'chainbreak',             'weight': '1.0'},
    ]

    as_names:List[str] = list()
    for res in structure.residues:
        if not res.is_ligand():
            continue

        chain_name:str = res.parent.name        
        rname:str=f"rs_{chain_name.lower()}"
        transform_name:str=f"dock_{chain_name.lower()}"
        elements.extend([
            {'parent': 'RESIDUE_SELECTORS', 'tag':'Index', 'name':rname, 'resnums':f"{res.idx}{res.parent.name}"},
            {'parent': 'RESIDUE_SELECTORS', 'tag':'CloseContact', 'name':f"{rname}_as", 'residue_selector':rname},
        ])

        as_names.append( f"{rname}_as" )

    elements.extend([    
        { 'parent': 'RESIDUE_SELECTORS', 'tag': 'Or', 'name':'as_selector', 'selectors':','.join(as_names)},
        { 'parent': 'RESIDUE_SELECTORS', 'tag': 'Not', 'name':'not_as_selector', 'selector':'as_selector'},
        {'parent':'MOVERS', 'tag':'FastRelax', 'name':'frelax', 'scorefxn':'hard_rep', 'cst_file':f"{Path(cst_file).absolute()}", 'ramp_down_constraints':'true', 'repeats':'5', 'child_nodes': [
                {'parent':'FastRelax', 'tag':'MoveMap', 'name':'full_enzyme', 'bb':'false', 'chi':'false', 'jump':'false','child_nodes':[
                        {'parent':'MoveMap','tag':'ResidueSelector', 'selector':'as_selector', 'bb':'true', 'chi':'true', 'bondangle':'true'},
                ]
        },]},
        {'parent':'PROTOCOLS', 'tag':'Add', 'mover_name':'frelax'},
    ])
    interface.rosetta.write_script(xml_file, elements)
#yapf: enable
    _parser = PDBParser()
    _parser.save_structure(pdb_file, structure)

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
        "-parser",
        f"   -protocol {Path(xml_file).absolute()}",
        "-out",
        f"   -file:scorefile 'score.sc'",
        "   -level 200",
        "   -nstruct 10",
        "   -overwrite",
        "   -path",
        f"       -all './rosetta_mm_min'",

    ])

    fs.write_lines(option_file, content)
    
    fs.safe_rmdir(f"{work_dir}/rosetta_mm_min/")
    fs.safe_mkdir(f"{work_dir}/rosetta_mm_min/")

    start_dir:str = Path(os.getcwd()).absolute()
    os.chdir( work_dir )
    interface.rosetta.run_rosetta_scripts([f"@{Path(option_file).name}"])
    os.chdir(start_dir)


    score_sc:str=f"{work_dir}/rosetta_mm_min/score.sc"
    assert Path(score_sc).exists()

    df: pd.DataFrame = interface.rosetta.parse_score_file(score_sc)
    df = df.sort_values(by='total_score').reset_index(drop=True)
    ref_stru = _parser.get_structure( f"{work_dir}/rosetta_mm_min/{df.iloc[0].description}.pdb" )
    stru_oper.update_residues(structure, ref_stru)

    for cst in constraints:
        cst.change_topology(structure)


#CURRENT    

