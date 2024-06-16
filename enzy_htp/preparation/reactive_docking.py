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
from enzy_htp._interface import (
    RosettaOptions,
    RosettaScriptsElement,
    RosettaScriptsProtocol
)
import enzy_htp.structure.structure_operation as stru_oper
from enzy_htp.structure.structure_constraint import (
    StructureConstraint, 
    CartesianFreeze
)

from enzy_htp.structure.structure_cluster import (
    cluster_structures,
    StructureCluster
)

import enzy_htp.chemical as chem
from enzy_htp.structure import PDBParser, Mol2Parser, Structure, Ligand, translate_structure, Atom
from enzy_htp.core import file_system as fs
from enzy_htp.quantum import single_point
from enzy_htp.quantum import optimize as qm_optimize 

def dock_reactants(structure: Structure,
                   ligands: List[Ligand],
                   constraints: List[StructureConstraint] = None,
                   n_struct: int = 100,
                   cst_energy: float = None,
                   use_qm: bool = True,
                   clash_cutoff: int = 3,
                   max_sasa_ratio:float = 0.60,
                   cluster_distance: float = 2.0,
                   contact_threshold:float=3.0,
                   box_size:float=50.0,
                   move_distance:float=40.0,
                   transform_angle:float=360,
                   transform_cycles:int=1000,
                   transform_repeats:int=3,
                   transform_temperature:int=5,
                   docking_fr_repeats:int=1,
                   min_fr_repeats:int=20,
                   grid_width:float=50.0,
                   rng_seed: int = 1996,
                   cluster_binding_dist:float=1.5,
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

    translate_structure(structure, end_naming='rosetta')
    
    interface.rosetta.parameterize_structure( structure, work_dir )

    for ligand in ligands:
        if id(ligand.root()) != id(structure):
            err_msg:str=f"The supplied ligand {ligand} is not a child of the supplied structure!"
            raise TypeError(err_msg)

    ligands.reverse()
    while ligands:
        ligand = ligands.pop()
        relevant_csts:List[StructureConstraint] = list()
        for cst in constraints:
            if cst.is_constraining(ligand):
                for ll in ligands:
                    if cst.is_constraining(ll):
                        break
                else:
                    relevant_csts.append(cst)

        opts = RosettaOptions()

        ### boiler plate options that are always the same
        opts['keep_input_protonation_state'] = True
        opts['auto_setup_metals'] = True
        opts['run:preserve_header'] = True
        opts["packing:ex1"] = True
        opts["packing:ex2"] = True
        opts["packing:ex2aro"] = True
        opts["packing:no_optH"] = False
        opts["packing:flip_HNQ"] = True
        opts["packing:ignore_ligand_chi"] = True
        opts['out:overwrite'] = True

        opts['run:constant_seed'] = True
        opts['run:jran'] = rng_seed
        opts['nstruct'] = n_struct
        opts['out:path:all'] = work_dir
        
        ### parameterization
        extra_res_fa:List[str] = list()
        stub_parent: str = os.path.expandvars(
            f"${config['rosetta.ROSETTA3']}/database/chemical/residue_type_sets/fa_standard/residue_types/protonation_states/")
        
        for stub in "GLU_P1 GLU_P2 LYS_D ASP_P1 TYR_D HIS_P ASP_P2".split():
            extra_res_fa.append(f"{stub_parent}/{stub}.params")

        if structure.data['rosetta_params']:
            extra_res_fa.extend( structure.data['rosetta_params'] )
            
        opts['extra_res_fa'] = ' '.join(map(lambda erf: f"'{erf}'", extra_res_fa))
        
        ### script variables
        opts.add_script_variable( 'ligand_chain', ligand.parent.name )
        opts.add_script_variable( 'ligand_idx', f"{ligand.idx}{ligand.parent.name}" )
        opts.add_script_variable( 'grid_width', grid_width )
        opts.add_script_variable( 'contact_threshold', contact_threshold )
        opts.add_script_variable( 'clash_cutoff', clash_cutoff )
        opts.add_script_variable( 'sasa_cutoff', max_sasa_ratio ) 
        opts.add_script_variable( 'cst_cutoff', cst_energy )
        opts.add_script_variable( 'box_size', box_size )
        opts.add_script_variable( 'move_distance', move_distance )
        opts.add_script_variable( 'transform_angle', transform_angle )
        opts.add_script_variable( 'transform_cycles', transform_cycles )
        opts.add_script_variable( 'transform_repeats', transform_repeats )
        opts.add_script_variable( 'transform_temperature', transform_temperature )
        opts.add_script_variable( 'fr_repeats', docking_fr_repeats )

        dock_ligand(structure,
                        ligand,
                        relevant_csts,
                        opts,
                        use_qm,
                        cluster_distance
                        )

    sp = PDBParser()

    opts.add_script_variable( 'fr_repeats', min_fr_repeats )
    opts.add_script_variable('ramp_constraints', False)
    opts['nstruct'] = 1 #Ooops
    mm_minimization(structure, constraints, opts)

    if use_qm:
        qm_minimization(structure, constraints, cluster_distance, False, work_dir)

    opts.add_script_variable('ramp_constraints', True)
    mm_minimization(structure, constraints, opts)

    if use_qm:
        qm_minimization(structure, [], cluster_distance, True, work_dir)

    translate_structure(structure, start_naming='rosetta')

    if not save_work_dir:
        _LOGGER.info(f"save_work_dir set to False! Deleting {work_dir}")
        fs.safe_rmdir( work_dir )


def mm_minimization(structure:Structure,
                constraints:List[StructureConstraint],
                opts:RosettaOptions,
                ) -> None:
    """
    """

    protocol = RosettaScriptsProtocol()
    protocol.add_residue_selector(
		'Index', name="ligand", resnums="%%ligand_idx%%"
    ).add_residue_selector(
		'CloseContact', name="ligand_active_site", residue_selector="ligand", contact_threshold="%%contact_threshold%%"
    ).add_residue_selector(
		'Not', name="not_ligand_active_site", selector="ligand_active_site"
    ).add_scorefunction(
        'ScoreFunction', name='hard_rep', weights='ligand'        
    ).add_mover(
        'FastRelax', name="frelax", scorefxn="hard_rep", cst_file="%%cst_file%%", repeats="%%fr_repeats%%", ramp_down_constraints="%%ramp_constraints%%", children=[
			('MoveMap', {'name':"full_enzyme", 'bb':"true", 'chi':"true", 'jump':"true", 'children':[
				('ResidueSelector', {'selector':"ligand_active_site",     'bb':"true", 'chi':"true", 'bondangle':"true" }),
				('ResidueSelector', {'selector':"not_ligand_active_site", 'bb':"false", 'chi':"true", 'bondangle':"false"})
                ]})]
    ).add_protocol(
        mover_name='frelax'
    )
    opts.add_script_variable('cst_file', 
        interface.rosetta.write_constraint_file( structure, constraints, opts['out:path:all'])
    )
    
    sele = list()
    for cst in constraints:
        for atom in cst.atoms:
            parent_residue = atom.parent
            sele.append(f"{parent_residue.idx}{parent_residue.parent.name}")

    opts.add_script_variable( 'ligand_idx',  ",".join(list(set(sele))))
    
    _LOGGER.info("Beginning Minimization step geometry sampling step...") #TODO(CJ):
    interface.rosetta.run_rosetta_scripts(
        structure,
        protocol,
        opts)

    df: pd.DataFrame = interface.rosetta.parse_score_file(opts['out:file:scorefile'], opts['out:path:all'])

    
    energy_key='total_score'

    infile:str=df.sort_values(by=energy_key).description.to_list()[0]
    
    _parser = PDBParser()
    ref_stru = _parser.get_structure(infile)
    stru_oper.update_residues(structure, ref_stru)

    for cst in constraints:
        cst.change_topology(structure)

def dock_ligand(structure:Structure,
                    ligand:Ligand,
                    constraints:List[StructureConstraint],
                    opts: RosettaOptions,
                    use_qm:bool,
                    cluster_distance:float=None
                    ) -> None:
    """TODO(CJ)""" 
    
    opts['qsar:grid_dir'] = str( Path(opts['out:path:all'] + "/qsar_grids/").absolute())
    fs.safe_mkdir(opts['qsar:grid_dir'])
    
    protocol = RosettaScriptsProtocol()
    protocol.add_scoring_grid(
        ligand_chain='%%ligand_chain%%', width='%%grid_width%%', name='grid',children=[
            ('ClassicGrid', {'grid_name':'classic', 'weight':'1.0'})]
    ).add_residue_selector(
        'Index', name='ligand', resnums='%%ligand_idx%%'
    ).add_residue_selector(
        'CloseContact', name='ligand_active_site', residue_selector='ligand', contact_threshold="%%contact_threshold%%"
    ).add_residue_selector(
        'Not', name="not_ligand_active_site", selector="ligand_active_site"
    ).add_scorefunction(
        'ScoreFunction', name='ligand_soft_rep', weights='ligand_soft_rep', children=[
			('Reweight', {'scoretype':"coordinate_constraint", 'weight':"1.0"}),
			('Reweight', {'scoretype':"atom_pair_constraint", 'weight':"1.0"}),
			('Reweight', {'scoretype':"angle_constraint", 'weight':"1.0"}),
			('Reweight', {'scoretype':"dihedral_constraint", 'weight':"1.0"}),
			('Reweight', {'scoretype':"chainbreak", 'weight':"1.0"})
        ]
    ).add_scorefunction(
        'ScoreFunction', name='hard_rep', weights='ligand'
    ).add_constraint_generator(
        'FileConstraintGenerator', name='add_cst', filename='%%cst_file%%'
    ).add_mover(
        'Transform', name='dock', chain='%%ligand_chain%%', box_size='%%box_size%%', 
            move_distance='%%move_distance%%',  angle="%%transform_angle%%", cycles="%%transform_cycles%%",
            repeats="%%transform_repeats%%", temperature="%%transform_temperature%%", grid_set="grid", 
            use_constraints="true", cst_fa_file="%%cst_file%%"
    ).add_mover(
        "ConstraintSetMover", name="add_csts", add_constraints="true", cst_file="%%cst_file%%" 
    ).add_mover(
        "ClearConstraintsMover", name="rm_csts"
    ).add_mover(
        'FastRelax', name="frelax", scorefxn="hard_rep", cst_file="%%cst_file%%", repeats="%%fr_repeats%%", children=[
            ('MoveMap', {'name':"full_enzyme", 'bb':"true", 'chi':"true", 'jump':"true", 'children':[
                ('ResidueSelector', {'selector':"ligand_active_site", 'bb':"true", 'chi':"true", 'bondangle':"true"}),
				('ResidueSelector', {'selector':"not_ligand_active_site", 'bb':"false", 'chi':"true", 'bondangle':"false"})
                ]}
            )]
    ).add_filter(
		'SimpleMetricFilter', name="clash_filter", comparison_type="lt_or_eq", cutoff="%%clash_cutoff%%", composite_action="any", children=[
			('PerResidueClashMetric', {'name':"clash", 'residue_selector':"ligand", 'residue_selector2':"ligand_active_site"})]
    ).add_filter(
        "SimpleMetricFilter", name="sasa_filter", comparison_type="lt_or_eq", cutoff="%%sasa_cutoff%%", composite_action="any", children=[
            ('SasaMetric', {'name':"sasa", 'residue_selector':"ligand"})  ]
    ).add_filter(
        'ConstraintScore', name="cst_filter", constraint_generators="add_cst", threshold="%%cst_cutoff%%" 
    ).add_protocol(
        mover_name="dock"
    ).add_protocol(
        filter="clash_filter"
    ).add_protocol(
        filter="sasa_filter"
    ).add_protocol(
        mover_name="frelax"    
    ).add_protocol(
        filter="cst_filter"    
    ).add_protocol(
        mover_name="rm_csts"
    )
    
    ligand_area:float=interface.pymol.get_ligand_area( structure, ligand )
    
    opts.add_script_variable('sasa_cutoff',  int(opts.get_script_variable('sasa_cutoff')*ligand_area))

    opts.add_script_variable('cst_file', 
        interface.rosetta.write_constraint_file( structure, constraints, opts['out:path:all'])
    )

    _LOGGER.info("Beginning RosettaLigand geometry sampling step...")
    interface.rosetta.run_rosetta_scripts(
        structure,
        protocol,
        opts)
       
    #TODO(CJ): add the clustering here

    df: pd.DataFrame = interface.rosetta.parse_score_file(opts['out:file:scorefile'], opts['out:path:all'])

    sfxn = None
    if use_qm:
        pass
        #TODO(CJ)
    
    _parser = PDBParser()
    structures:List[Structure] = list(map(lambda dd: _parser.get_structure( dd ), df.description ) )
    
    clusters:List[StructureCluster] = cluster_structures( structures, 'polymer.protein', f"resn {ligand.name}", 1.0 ) #TODO(CJ): update params
    for cc in clusters:
        interface.rosetta.score( cc, opts )

    ref_stru = sorted(
        clusters,
        key=lambda clust: clust.average_score()
    )[0].lowest_energy_structure()

    stru_oper.update_residues(structure, ref_stru)

    for cst in constraints:
        cst.change_topology(structure)

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
                freeze_ligands:bool,
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
   
    to_freeze:List[Atom] = list()

    for res in structure.residues:
        if res.is_canonical():
            continue

        for atom in res.atoms:
            atom.charge = 0.0

            if freeze_ligands and atom.element != 'H':
                to_freeze.append(atom)
    
    translate_structure(structure, start_naming='rosetta')
    es = qm_optimize(structure,
            engine="xtb",
            constraints=constraints + [CartesianFreeze(structure.backbone_atoms() + to_freeze )],
            regions=[as_sele],
            region_methods=[chem.QMLevelOfTheory(basis_set='',method='GFN2', solvent='water', solv_method='ALPB')],
            parallel_method=None)[0]

    translate_structure(structure, end_naming='rosetta')

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

    _LOGGER.info(f"Beginning qm energy evaluation.")
    as_sele:str = get_active_site_sele(structure, cluster_cutoff)
    qm_energy = []

    _parser = PDBParser()
    for i, row in df.iterrows():

        energy:float = None
        _df_stru = _parser.get_structure( row.description )
        translate_structure(_df_stru, start_naming='rosetta')
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
