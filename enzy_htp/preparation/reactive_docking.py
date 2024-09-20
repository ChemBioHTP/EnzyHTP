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
                   stereo_constraints: List[StructureConstraint] = None,
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

    if stereo_constraints is None:
        stereo_constraints = []

    fs.safe_mkdir(work_dir)
    fs.safe_mkdir("./snapshots")
    
    sp = PDBParser()

    sp.save_structure("./snapshots/initial_structure.pdb", structure)

    translate_structure(structure, end_naming='rosetta')
    
    interface.rosetta.parameterize_structure( structure, work_dir )

    for ligand in ligands:
        if id(ligand.root()) != id(structure):
            err_msg:str=f"The supplied ligand {ligand} is not a child of the supplied structure!"
            raise TypeError(err_msg)

    ligands.reverse()
    lrk = set() 
    while ligands:
        ligand = ligands.pop()
        relevant_csts:List[StructureConstraint] = list()
        for cst in constraints + stereo_constraints:
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
        stub_parent: str = os.path.expandvars(
            f"${config['rosetta.ROSETTA3']}/database/chemical/residue_type_sets/fa_standard/residue_types/protonation_states/")
        
        for stub in "GLU_P1 GLU_P2 LYS_D ASP_P1 TYR_D HIS_P ASP_P2".split():
            opts.add_extra_res_fa(  f"{stub_parent}/{stub}.params" )

        if structure.data['rosetta_params']:
            for pp in structure.data['rosetta_params']:
                opts.add_extra_res_fa( pp )
            
        #opts['extra_res_fa'] = ' '.join(map(lambda erf: f"'{erf}'", extra_res_fa))
        
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
        
        lrk.add( ligand.key() )
        dock_ligand(structure,
                        ligand,
                        relevant_csts,
                        opts,
                        use_qm,
                        cluster_distance
                        )
    

    sp.save_structure("./snapshots/docked_structure.pdb", structure)
    
    opts.add_script_variable( 'fr_repeats', min_fr_repeats )
    opts.add_script_variable('ramp_constraints', False)
    opts['nstruct'] = 1 #Ooops
    
    mm_minimization(structure, constraints, opts, True)
    sp.save_structure("./snapshots/mm_structure_01.pdb", structure)

    if use_qm:
        #TODO(CJ): this will be the xtb only optimization
        qm_minimization_xtb(structure, constraints, cluster_distance, False, work_dir, lrk)
        sp.save_structure("./snapshots/qm_structure_01.pdb", structure)

    opts.add_script_variable('ramp_constraints', True)
    mm_minimization(structure, constraints, opts, False)
    sp.save_structure("./snapshots/mm_structure_02.pdb", structure)

    if use_qm:
        qm_minimization(structure, constraints, cluster_distance, True, work_dir, lrk)
        sp.save_structure("./snapshots/qm_structure_02.pdb", structure)

    translate_structure(structure, start_naming='rosetta')

    if not save_work_dir:
        _LOGGER.info(f"save_work_dir set to False! Deleting {work_dir}")
        fs.safe_rmdir( work_dir )


def mm_minimization(structure:Structure,
                constraints:List[StructureConstraint],
                opts:RosettaOptions,
                bb_flex:bool
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
			('MoveMap', {'name':"full_enzyme", 'bb':"true" if bb_flex else "false", 'chi':"true", 'jump':"true", 'children':[
				('ResidueSelector', {'selector':"ligand_active_site",     'bb': "true" if bb_flex else 'false', 'chi':"true", 'bondangle':"true" }),
				('ResidueSelector', {'selector':"not_ligand_active_site", 'bb': "true" if bb_flex else "false", 'chi':"true", 'bondangle':"false"})
                ]})]
    ).add_protocol(
        mover_name='frelax'
    )
    opts.add_script_variable('cst_file', 
        interface.rosetta.write_constraint_file( structure, constraints, work_dir=opts['out:path:all'])
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
        opts,
        prefix="mm_minimization"
        )

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
        filter="cst_filter"    
    ).add_protocol(
        mover_name="rm_csts"
    ).add_protocol(
        mover_name="frelax"    
    )
    
    ligand_area:float=interface.pymol.get_ligand_area( structure, ligand )
    
    opts.add_script_variable('sasa_cutoff',  int(opts.get_script_variable('sasa_cutoff')*ligand_area))

    opts.add_script_variable('cst_file', 
        interface.rosetta.write_constraint_file( structure, constraints, work_dir=opts['out:path:all'])
    )

    _LOGGER.info("Beginning RosettaLigand geometry sampling step...")
    interface.rosetta.run_rosetta_scripts(
        structure,
        protocol,
        opts,
        prefix="docking"
        )

    df: pd.DataFrame = interface.rosetta.parse_score_file(opts['out:file:scorefile'], opts['out:path:all'])

    _parser = PDBParser()
    structures:List[Structure] = list(map(lambda dd: _parser.get_structure( dd ), df.description ) )
    
    clusters:List[StructureCluster] = cluster_structures( structures, 'polymer.protein', f"resn {ligand.name}", 1.0 ) #TODO(CJ): update params
    
    score_clusters(
        clusters,
        structure,
        cluster_distance,#TODO(CJ): I think this comes from the options?
        opts,
        use_qm
    )

    ref_stru = sorted(
        clusters,
        key=lambda clust: clust.average_score()
    )[0].lowest_energy_structure()

    stru_oper.update_residues(structure, ref_stru)

    for cst in constraints:
        cst.change_topology(structure)

def get_active_site_sele(structure: Structure, distance_cutoff:float, fmt:str='pymol', bridge_gaps:bool=False, constraints=None, lrk=None) -> str:
    """Creates a pymol-compatible sele for the active site of the supplied Structure(). Basic structure
    is to select all residues within the specified cutoff of the Ligand()'s in the structure. Metal ions up to
    2*distance_cutoff from the Ligand()'s are also selected.
    
    Args:
        structure: The Structure() to be used as a template for the active site.
        distance_cutoff: The cutoff in Angstroms for a Residue() to be included in the active site.  
        
    Returns:
        Selection string in pymol format which defines the enzyme active site.        
    """
    #TODO(CJ): make this work with constraints
    _LOGGER.info("Analyzing enzyme active site...")
    #for res in structure.residues:
    ligand_residue_keys = set()
    if lrk is None:
        for res in structure.residues:
            if res.is_ligand():
                ligand_residue_keys.add(res.key())
    else:
        ligand_residue_keys = lrk

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
        result.add((row['chain'], int(row['resi'])))

    if constraints is not None:
        for cst in constraints:
            for atom in cst.atoms:
                result.add( atom.parent.key() )

    if bridge_gaps:
        result = list(result)
        result.sort()
        bridge = list()
        for r1, r2 in zip(result[:-1], result[1:]):
            if r1[0] != r2[0]:
                continue
            
            if abs(int(r1[1]) -  int(r2[1])) == 2:
                _LOGGER.info(f"Bridging gap between {r1[0]}.{r1[1]} and {r2[0]}.{r2[1]}")
                bridge.append((r1[0], int(r1[1])+1))
            
            result.extend( bridge )

    _LOGGER.info(f"Found {len(result)} residues within {distance_cutoff} angstroms of reactants!")
    
    if fmt == 'pymol':
        return " or ".join(map(
            lambda rr: f"( chain {rr[0]} and resi {rr[1]})",
            result
        ))
    elif fmt == 'rosetta':
        return ",".join(map(
            lambda rr: f"{rr[1]}{rr[0]}", result
        ))
    else:
        assert False

def qm_minimization_xtb(structure:Structure,
                constraints:List[StructureConstraint],
                cluster_distance:float,
                freeze_ligands:bool,
                work_dir:str,
                lrk=None
                ) -> None:


    as_sele:str = get_active_site_sele(
                structure,
                cluster_distance,
                'pymol',
                bridge_gaps=True,
                constraints=constraints,
                lrk=lrk)

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
                engine='xtb',
                constraints=constraints + [CartesianFreeze(structure.backbone_atoms() + to_freeze)],
                regions=[as_sele],
                region_methods=[chem.QMLevelOfTheory(basis_set='', method='GFN2', solvent='water', solv_method='ALPB')],
                parallel_method=None)[0]

    translate_structure(structure, end_naming='rosetta')

def qm_minimization(structure:Structure,
                constraints:List[StructureConstraint],
                cluster_distance:float,
                freeze_ligands:bool,
                work_dir:str,
                lrk=None) -> None:
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
    #TODO(CJ): add to freeze stuff for ligands
    #TODO(CJ): update to the new version with better capping and CB freezing

    protocol = RosettaScriptsProtocol()
    options = RosettaOptions()
    protocol.add_residue_selector(
        "Index", name="active_site", resnums=get_active_site_sele(structure, cluster_distance, 'rosetta', bridge_gaps=True, constraints=constraints, lrk=lrk)
    ).add_scorefunction(
        "ScoreFunction", name="r15", weights="ref2015"
    ).add_scorefunction(
        "ScoreFunction", name="qm_region", children=[
            ( 'Reweight',  {'scoretype': 'orca_qm_energy', 'weight':'1.0'}),
            ( 'Set', {'orca_processes':'8'}),
            ( 'Set', {'orca_path':config['rosetta.ORCA_DIR'] }),
            ( 'Set', {'orca_memory_megabytes':"3000" }),
            ( 'Set', {'orca_electron_correlation_treatment':"XTB" }),
            ( 'Set', {'rosetta_orca_bridge_temp_directory':"xtb_temp" }),
            ( 'Set', {'orca_geo_opt_max_steps':"1500" }),
            ( 'Set', {'clean_rosetta_orca_bridge_temp_directory':"false" }),
            ( 'Set', {'orca_deduce_charge':"true" }),
        ]
    ).add_scorefunction(
        'MultiScoreFunction', name="combo_sfxn", children=[
            ('SimpleCombinationRule', {}),
            ('Region', {'scorefxn':'qm_region', 'residue_selector':'active_site', 'children':[
                ('CappedBondResolutionRule', {})
            ]}),
            ('Region', {'scorefxn':'r15', 'children':[
                ('SimpleBondResolutionRule', {})
            ]}),
        ]
    ).add_mover(
        'OrcaQMGeometryOptimizationMover', name='qm_opt',
            freeze_backbone_atoms='true',
            msfxn_name="combo_sfxn",
            msfxn_freeze_noncommon_atoms="true",
            clean_rosetta_orca_bridge_temp_directory='false',
            geo_opt_max_steps="1500", deduce_charge="true",
            orca_memory_megabytes="20000",
            immobilize_h_bond_lengths="true",
            optimization_convergence="LOOSEOPT" #TODO(CJ)
    ).add_mover(
        'ConstraintSetMover', name="add_cst",
            add_constraints="true",
            cst_file="%%cst_file%%"
    ).add_protocol(
        mover_name='add_cst'
    ).add_protocol(
        mover_name='qm_opt'
    )
    
    options['overwrite'] = True
    
    interface.rosetta.parameterize_structure( structure, work_dir )

    for prm in structure.data['rosetta_params']:
        options.add_extra_res_fa( prm )
    
    options.add_script_variable('cst_file', interface.rosetta.write_constraint_file(structure, constraints, 'ORCA_FROZEN', work_dir=work_dir))
    
    interface.rosetta.run_rosetta_scripts( structure, protocol, options)
    
    df: pd.DataFrame = interface.rosetta.parse_score_file(options['out:file:scorefile'], options['out:path:all'])
    
    infile:str=df.iloc[0].description
    
    _parser = PDBParser()
    
    ref_stru = _parser.get_structure(infile)
    stru_oper.update_residues(structure, ref_stru)
    
    for cst in constraints:
        cst.change_topology( structure )


def score_clusters(
    clusters,
    structure,
    cluster_cutoff,
    opts,
    use_qm
    ) -> None:
    #TODO(CJ)
    local_opts = deepcopy( opts )
    local_opts['nstruct'] = 1
    for cc in clusters:
        sfxn=None
        protocol=None

        if use_qm:
            protocol = RosettaScriptsProtocol()
            protocol.add_scorefunction(
                "ScoreFunction", name='total_region', weights='ref2015'
            ).add_scorefunction(
                "ScoreFunction", name="qm_region", children=[
                    ('Reweight', {'scoretype':"orca_qm_energy", 'weight':"1.0" }),
                    ('Set', {'orca_path':config['rosetta.ORCA_DIR'] }),
                    ('Set', {'orca_processes':"8" }),
                    ('Set', {'orca_memory_megabytes':"3000" }),
                    ('Set', {'rosetta_orca_bridge_temp_directory':"xtb_temp" }),
                    ('Set', {'orca_electron_correlation_treatment':"XTB" }),
                    ('Set', {'orca_geo_opt_max_steps':"1500" }),
                    ('Set', {'clean_rosetta_orca_bridge_temp_directory':"true" }),
                    ('Set', {'orca_deduce_charge':"true" })
                ]
            ).add_residue_selector(
                'Index', name='ligand', resnums=get_active_site_sele(structure, cluster_cutoff, 'rosetta')
            )
            sfxn = RosettaScriptsElement(
                'MultiScoreFunction', name='sfxn', children=[
                    ('SimpleCombinationRule',{}),
                    ('Region', {'scorefxn':'qm_region', 'residue_selector': 'ligand', 'children':[
                        ('CappedBondResolutionRule',)
                    ]}),
                    ('Region', {'scorefxn':'total_region', 'children':[
                        ('SimpleBondResolutionRule',)
                    ]})
                ]
            )

        interface.rosetta.score( 
            cc, 
            opts=local_opts,
            protocol=protocol,
            score_fxn=sfxn,
            prefix="xtb_score"
        )

