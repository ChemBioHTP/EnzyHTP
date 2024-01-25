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
from enzy_htp.structure.structure_constraint import StructureConstraint, CartesianFreeze

import enzy_htp.chemical as chem
from enzy_htp import mutation as mm
from enzy_htp.structure import PDBParser, Mol2Parser, Structure, Ligand
from enzy_htp.core import file_system as fs
from enzy_htp.quantum import single_point
from enzy_htp.quantum import optimize as qm_optimize 

def dock_reactants(structure: Structure,
                   constraints: List[StructureConstraint] = None,
                   n_struct: int = 100,
                   cst_energy: float = None,
                   use_qm: bool = True,
                   freeze_alphafill:bool = True,
                   clash_cutoff: int = 3,
                   max_sasa_ratio:float = 0.60,
                   cluster_distance: float = 2.0,
                   rng_seed: int = 1996,
                   work_dir: str = None,
                   save_work_dir: bool = True,
                   ) -> None:
    """Takes a Structure() containing Ligand() objects and tries to create a complex with an optimized geometry 
    consistent with that described in the supplied constraints. Does all work inplace on the supplied Structure().
    Below is the workflow of the function:

    1. Generate constrained geometries with RosettaLigand()
    2. Filter down and select a target geometry
        a. Filter out if too many clashes
        b. Filter out if ligands have too much SASA
        c. Filter out if ligands do not satisfy constraints
        d. Select geometry with:
            i. Lowest QM energy, if use_qm is True
            ii. Lowest RosettaLigand energy, if use_qm is False
    3. Rosetta FastRelax minimization 1
        a. use constraints
        b. allow backbone flexibility in active site
    4. constrained QM minimization of active site with xtb (if use_qm is True)
    5. Rosetta FastRelax minimization 2
        a. ramp down constraints
        b. oly sidechain flexibility in active site
    6. un-constrained QM minimization of active site with xtb (if use_qm is True)
    
    Args:
        structure: The Structure() to perform reactive docking on. Contains both Ligand()'s and protein chains.
        constraints: A List[StructureConstraint] describing the reaction geometry.
        n_struct: How many geometries should we make with RosettaLigand? Default is 100.
        cst_energy: The energy penalty cutoff for keeping a given geometry. In REU and assuming Rosetta constraint scoring.
        use_qm: Shoud QM be used for both geometry filtration and geometry minimization? Default is True.
        freeze_alphafill: If a Ligand() was placed with AlphaFill, should it be frozen during initial geometry sampling. Default is True.
        clash_cutoff: How many heavy atom clashes are allowed in a given geometry? Default is 3.
        max_sasa_ratio: What percentage of a Ligand()'s surface area can be solvent accessible?  Default is 0.60.
        cluster_distance: How close should a residue be to the Ligand()'s to be included in the QM region? Default is 2.0 Angstroms.
        rng_seed: Random number generation integer seed. Default value is 1996.
        work_dir: Where is the work being done/where should the scratch files be made?
        save_work_dir: Should temporary files be saved? Default is True.

    Returns:
        Nothing.
    """

    if work_dir is None:
        work_dir = config["system.SCRATCH_DIR"]

    if cst_energy is None:
        cst_energy = len(constraints)*2000.0

    fs.safe_mkdir(work_dir)

    interface.rosetta.rename_atoms(structure)
    
    param_files:List[str] = generate_ligand_params(structure, work_dir)
    
    geometry_df:pd.DataFrame = generate_geometries(structure, constraints, n_struct, clash_cutoff, param_files, freeze_alphafill, rng_seed, work_dir )

    select_geometry(structure, constraints, geometry_df, max_sasa_ratio, cst_energy, cluster_distance, use_qm)

    if use_qm:
        qm_minimization(structure, constraints, cluster_distance, work_dir)

#    mm_minimization(structure, cluster_distance, constraints, param_files, rng_seed, work_dir)
#
#    if use_qm:
#        qm_minimization(structure, [], cluster_distance, work_dir)

    if not save_work_dir:
        _LOGGER.info(f"save_work_dir set to False! Deleting {work_dir}")
        fs.safe_rmdir( work_dir )


def create_docking_xml(stru:Structure, cst_file:str, freeze_alphafill:bool, clash_cutoff:int, work_dir: str) -> str:
    """Creates the docking .xml RosettaScripts File that will be used for creating the reactive complex geometries.

    Args:
        stru: The Structure() that the geometry sampling will occur on.
        cst_file: The .cst file containing the atom constraints which define the system.
        freeze_alphafill: Should Ligand()'s that were placed by AlphaFill be frozen?
        clash_cutoff: How many heavy atom clashes are allowed for each geometry?
        work_dir: Where the temporary file will be written.

    Returns:
        The relative filepath of the RosettaScripts .xml file.
    """
    fname: str = f"{work_dir}/docking_script.xml"
    
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
    grid_index:int = 1
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
        if freeze_alphafill and res.placement_method == 'alphafill':
            continue
        
        grid_name:str=f'grid_{grid_index}'
        grid_index += 1
        dock:str = f"dock_{chain_name}"
        clash_metric:str=f"clash_{chain_name.lower()}"
        clash_filter:str=f"{clash_metric}_filter"
        elements.extend([{'parent': 'ROSETTASCRIPTS', 'tag': 'SCORINGGRIDS', 'ligand_chain': chain_name, 'width': '20', 'append_elements_only': True, 'name': grid_name, 'child_nodes': [
                    {
                        'parent': 'SCORINGGRIDS',
                        'tag': 'ClassicGrid',
                        'grid_name': 'classic',
                        'weight': '1.0'
                    },
                ]
            },
            { 'parent': 'MOVERS', 'tag': 'Transform', 'name': dock, 'chain': chain_name, 'box_size': str(50), 'move_distance': '20.0', 'angle': '360', 'cycles': '1000', 'repeats': '3', 'temperature': '5',
                'grid_set': grid_name, 'use_constraints':'true', 'cst_fa_file':str(Path(cst_file).absolute())},
            {'parent': 'PROTOCOLS', 'tag': 'Add', 'mover_name':  dock },
            {'parent': 'FILTERS', 'tag':'SimpleMetricFilter', 'name':clash_filter, 'comparison_type':'lt_or_eq', 'cutoff':f"{clash_cutoff}", 'composite_action':'any', 'child_nodes':[
                {'parent':'SimpleMetricFilter', 'tag':'PerResidueClashMetric', 'name':clash_metric, 'residue_selector':rname, 'residue_selector2': f"{rname}_as"}
            ]},
            {'parent':'PROTOCOLS', 'tag':'Add', 'filter':clash_filter}
            ])

    elements.extend([
        { 'parent': 'RESIDUE_SELECTORS', 'tag': 'Or', 'name':'as_selector', 'selectors':','.join(as_names)},
        { 'parent': 'RESIDUE_SELECTORS', 'tag': 'Not', 'name':'not_as_selector', 'selector':'as_selector'},
        { 'parent': 'MOVERS', 'tag': 'FastRelax', 'name': 'frelax', 'scorefxn': 'hard_rep', 'cst_file':f"{Path(cst_file).absolute()}", 'child_nodes': [
                {'parent':'FastRelax', 'tag':'MoveMap', 'name':'full_enzyme', 'bb':'true', 'chi':'true', 'jump':'false','child_nodes':[
                        {'parent':'MoveMap','tag':'ResidueSelector', 'selector':'as_selector', 'bb':'true', 'chi':'true', 'bondangle':'true'},
                        {'parent':'MoveMap', 'tag':'ResidueSelector', 'selector':'not_as_selector', 'bb':'false', 'chi':'true', 'bondangle':'false'},]}]
        },
        {'parent':'PROTOCOLS', 'tag':'Add', 'mover_name':'frelax'}
    ])

    interface.rosetta.write_script(fname, elements)
    _LOGGER.info(f"Saved new RosettaScripts .xml file at {fpath.absolute()}!")
    return fname


def create_docking_options_file(pdb_file: str,
                      xml_file: str,
                      param_files: List[str],
                      work_dir: str,
                      rng_seed: int,
                      n_struct: int,) -> str:
    """Makes the docking_options.txt file that the docking run will actually use. This function DOES NOT make any checks to the inputs.
    
    Args:
        pdb_file: The .pdb file (with constraints) to use. 
        xml_file: The validated RosettaScripts .xml file to be used for docking.
        param_files: The list() of reactant .params files.
        work_dir: The working directory 
        rng_seed: rng seed to be used during dcking.
        n_struct: Number of strutures to make as an int().

    Returns:
        Path to the docking_options.txt file with all the Rosetta options.
    """
    _LOGGER.info("Beginning creation of options file for Rosetta...")
    content: List[str] = [
        "-keep_input_protonation_state", 
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

    _LOGGER.info(f"Wrote the below settings to {fname}:")
    for ll in content:
        _LOGGER.info(f"\t{ll}")
    fs.write_lines(fname, content)

    option_file = fname.absolute()

    return str(fname)


def generate_geometries(structure: Structure,
            constraints: List[StructureConstraint],
            n_struct:int,
            clash_cutoff:int,
            param_files:List[str],
            freeze_alphafill:bool,
            rng_seed:int,
            work_dir:str,
            ) -> pd.DataFrame:
    """Geometry generation engine that leverages RosettaLigand to create the specified number of geometries. Also filters out
    geometries with too many clashes. Returns the contents of the Rosetta score.sc file as a pandas DataFrame where description
    is subbed out to be the absolute path of the corresponding .pdb file. 

    Args:
        structure: The Structure() object to generate geometries with.
        constraints: The List[StructureConstraint] that defines the reactive complex geometry.
        n_struct: How many geometries should be created?
        clash_cutoff: How many clashes are allowed to be in each geometry before they are removed.
        param_files: A List[str] containing .params files for use in Rosetta.
        freeze_alphafill: Should Ligand()'s placed with AlphaFill be frozen during geometry generation.
        rng_seed: The integer random number generation seed to use during geometry generation.
        work_dir: The path to where all the work should be done.

    Returns:
        A pandas DataFrame that contains the contents of the score.sc file. The description column contains 

    """
    start_pdb:str = f"{work_dir}/start.pdb"
    parser = PDBParser()
    parser.save_structure(start_pdb, structure)

    cst_file:str = interface.rosetta.write_constraint_file(structure, constraints, work_dir) #TODO(CJ): look at this; wrong constraint types!!

    xml_file:str = create_docking_xml(structure, cst_file, freeze_alphafill, clash_cutoff, work_dir)  #TODO(CJ): going to overhaul this

    options_file: str = create_docking_options_file(start_pdb,
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

    scores_file:str = str((opt_path.parent / "complexes/score.sc").absolute())
    
    fs.check_file_exists(scores_file, exit_script = False )

    df: pd.DataFrame = interface.rosetta.parse_score_file(scores_file)

    df['description'] = df.apply(lambda row: f"{opt_path.parent}/complexes/{row.description}.pdb", axis=1)

    df['selected'] = True

    _LOGGER.info("Completed RosettaLigand geometry sampling!")

    return df

def evaluate_geometry_SASA(df: pd.DataFrame, max_sasa_ratio:float) -> None:
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
        raise TypeError()

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


def evaluate_geometry_csts(df: pd.DataFrame, csts: List[StructureConstraint], cst_cutoff: int) -> None:
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
    

def generate_ligand_params(stru:Structure, work_dir:str) -> List[str]:
    """Given the input Structure(), parameterize everything needed to use the Ligand()'s in Rosetta.
    
    Args:
        stru: The Structure() to parameterize.
        work_dir: Where temporary files will be saved.

    Returns:
        A List[str] with ligand .params files.         
    """
    #TODO(CJ): probably move this to the RosettaInterface
    _LOGGER.info("Beginning preparation of each reactant...")
    param_files: List[str] = list()
   
    parser = Mol2Parser()
    for res in stru.residues:
        conformers:str=None
        if type(res) != Ligand:
            continue
        
        _LOGGER.info(f"Detected residue {res.name} in chain {res.parent.name}...")

        if res.net_charge is None:
            res.net_charge = interface.bcl.calculate_formal_charge( res )

        param_file:str = interface.rosetta.parameterize_ligand(res, charge=res.net_charge, work_dir=work_dir)
        
        _LOGGER.info(f"Information for reactant {res.name}:")
        _LOGGER.info(f"\tparam file: {param_file}")
        _LOGGER.info(f"\tcharge: {res.net_charge}")
        _LOGGER.info(f"\tnum conformers: {res.n_conformers()}")
        
        param_files.append(param_file)
    _LOGGER.info("Finished reactant preparation!")

    return param_files


def get_active_site_sele(structure: Structure, distance_cutoff:float) -> str:
    """Creates a pymol-compatible sele for the active site of the supplied Structure(). Basic structure
    is to select all residues within the specified cutoff of the Ligand()'s in the structure. Metal ions up to
    2*distance_cutoff from the Ligand()'s are also selected.
    
    Args:
        structure: The Structure() to be used as a template for the active site.
        distance_cutoff: The cutoff in Angstroms for a Residue() to be included in the active site.  
        
    Returns:
        Selection string in pymol format which defines the enzyme active site.        
    """

    _LOGGER.info("Analyzing enzyme active site...")
    #for res in structure.residues:
    ligand_residue_keys = set()
    for res in structure.residues:
        if res.is_ligand():
            ligand_residue_keys.add(res.key())

    parser = PDBParser()
    start_pdb:str=f"{config['system.SCRATCH_DIR']}/active_site_selection.pdb"
    parser.save_structure(start_pdb, structure)
   
    ligand_sele:str = f"byres ( {'or '.join(map(lambda lrk: f'(all within {distance_cutoff} of chain {lrk[0]} and resi {lrk[1]})', ligand_residue_keys))}  )"
    ligand_sele:str = " or ".join(map(
        lambda lrk: f"((byres all within {distance_cutoff:.2f} of (chain {lrk[0]} and resi {lrk[1]})) or (metals within  {2*distance_cutoff:.2f} of( chain {lrk[0]} and resi {lrk[1]})))",
        ligand_residue_keys
    ))

    session = interface.pymol.new_session()
    df = interface.pymol.collect(session, start_pdb, "chain resi".split(), sele=ligand_sele)

    fs.safe_rm(start_pdb)

    result = set()
    for i, row in df.iterrows():
        result.add((row['chain'], row['resi']))
    
    _LOGGER.info(f"Found {len(result)} residues within {distance_cutoff} angstroms of reactants!")
    
    return " or ".join(map(
        lambda rr: f"( chain {rr[0]} and resi {rr[1]})",
        result
    ))

def qm_minimization(structure:Structure,
                constraints:List[StructureConstraint],
                cluster_distance:float,
                work_dir:str) -> None:
    """Performs QM minimization of the enzyme active site using xtb. Assumes that backbone atoms of the Residue()'s should be 
    frozen. Is capable of converting supplied constraints to xtb format. Updates coordinates in place. 
    
    Args:
        structure: The Structure() object to optimize.
        constraints: The List[StructureConstraint] which define the enzyme active site geometry.
        cluster_distance: The cutoff in angstroms for inclusion of a Residue() in the active site.
        work_dir: The temporary directory to do all work.

    Returns:
        Nothing.
    """

    as_sele = get_active_site_sele(structure, cluster_distance)
    
    for res in structure.residues:
        if res.is_canonical():
            continue
        
        for atom in res.atoms:
            atom.charge = 0.0
   
    es = qm_optimize(structure,
            engine="xtb",
            constraints=constraints + [CartesianFreeze(structure.backbone_atoms())],
            regions=[as_sele],
            region_methods=[chem.QMLevelOfTheory(basis_set='',method='GFN2', solvent='water', solv_method='ALPB')],
            parallel_method=None)[0]

def evaluate_geometry_qm_energy(df: pd.DataFrame, structure: Structure, cluster_cutoff: float) -> None:
    """Aids in ranking and selection of candidate geometries through a semi-empirical QM single point energy
    calculation with xtb. Creates a capped active site by using a specified cluster_cutoff parameter to specify
    the enzyme's active site.

    Args:
        df: The geometry DataFrame containing all information  
        structure: The reference Structure() in use.
        cluster_cutoff: The cutoff in Angstroms for a Residue() to be included in the QM region. 

    Returns:
        Nothing.        
    """
    if not df.selected.sum():
        _LOGGER.error("No geometries are still selected!")
        raise TypeError()

    _LOGGER.info(f"Beginning qm energy evaluation. {df.selected.sum()} geometries still selected...")
    as_sele:str = get_active_site_sele(structure, cluster_cutoff)
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
        es = single_point(
            _df_stru,
            engine='xtb',
            region_methods=[chem.QMLevelOfTheory(basis_set='',method='GFN2', solvent='water', solv_method='ALPB')],
            parallel_method=None,
            regions=[as_sele])
        qm_energy.append(es[0].energy_0)

    df['qm_energy'] = qm_energy
    _LOGGER.info("Finished qm energy evaluation!")

def select_geometry( structure:Structure,
                        constraints:List[StructureConstraint],
                        df:pd.DataFrame,
                        max_sasa_ratio:float,
                        cst_energy:float,
                        cluster_distance:float,
                        use_qm:bool) -> None:
    """Method that takes the results from the generate_geometries() function and selects the geometry to move forward with. The coordinates
    from the selected geometry are applied directly to the supplied Structure(), in place. The function filters out candidate geometries
    if the ligands have too much SASA (as specified by the user) or the actual geometry deviates too much from the geometry specified by the
    StructureConstraint. The selected geometry is determined by an energy ranking. If use_qm=True, then an active site calculation is performed
    with xtb, else the total energy metric from the RosettaLigand run is used.

    Args:
        structure: The Structure() to apply the selected geometry to.
        constraints: A List[StructureConstraint] which define the reacting geometry.
        df: The pd.DataFrame from the generate_geometries() function.
        max_sasa_ratio: The maximum SASA ratio the ligands may have. 
        cst_energy: The maximum constraint energy penalty allowed in the geometry.
        cluster_distance: How far away can Residue()'s be in the active site?
        use_qm: Should QM energies be used to finally select the geometry? 

    Returns:
        Nothing.
    """
    evaluate_geometry_SASA(df, max_sasa_ratio)#TODO(CJ): put this back in 

    evaluate_geometry_csts(df, constraints, cst_energy)
   
    energy_key:str=None
    
    if use_qm:
        evaluate_geometry_qm_energy(df, structure, cluster_distance)
        energy_key = 'qm_energy'
    else:
        energy_key = 'total_score'
   
    if not df.selected.sum():
        _LOGGER.error("No geometries satisfy all selection criteria!")
        raise TypeError()

    infile:str=df[df.selected].sort_values(by=energy_key).description.to_list()[0]

    _parser = PDBParser()
    ref_stru = _parser.get_structure(infile)
    stru_oper.update_residues(structure, ref_stru)

    for cst in constraints:
        cst.change_topology(structure)

def mm_minimization(structure:Structure, cluster_distance:float, constraints:List[StructureConstraint], param_files:List[str], rng_seed:int, work_dir:str) -> None:
    """Uses Rosetta's FastRelax protocol to perform a minimization with constraints that ramp-down over the course of the minimization. Updates the coordinates
    of the Structre() inplace. Residues are given differing levels of flexibility where there is no backbone or sidechain for all of the enzyme except for the
    active site. In the active site there is sidechain and bond angle flexibility. The geometry with the lowest total_score is selected.

    Args:
        structure: The Structure() to operate on.
        cluster_distance: The distance cutoff in Angstroms to be an "active site" residue.
        constraints: The List[StructureConstraint] which define the reactive complex geometry.
        param_files: A List[str] of .params files for use in Rosetta.
        rng_seed: An int() seed for Rosetta's random number generator.
        work_dir: The directory where temporary files and directories will be written. 

    Returns:
        Nothing.
    """

    cst_file:str = interface.rosetta.write_constraint_file(structure, constraints, work_dir) 
    xml_file:str = f"{work_dir}/mm_script.xml"
    pdb_file:str = f"{work_dir}/temp_mm_min.pdb"
    option_file:str = f"{work_dir}/mm_min_options.txt"
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
                        {'parent':'MoveMap','tag':'ResidueSelector', 'selector':'as_selector', 'bb':'false', 'chi':'true', 'bondangle':'true'},
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

