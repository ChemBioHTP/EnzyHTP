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
    RosettaScriptsProtocol,
    RosettaScriptsEngine
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
from enzy_htp.core import (
    file_system as fs,
    job_manager as jm
)
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
                   min_fr_repeats:int=10,
                   grid_width:float=50.0,
                   rng_seed: int = 1996,
                   cluster_binding_dist:float=1.5,
                   work_dir: str = None,
                   save_work_dir: bool = True,
                   cpu_config = None,
                   qm_config = None
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
    
    wstru = deepcopy( structure )

    for ligand in ligands:
        if id(ligand.root()) != id(structure):
            err_msg:str=f"The supplied ligand {ligand} is not a child of the supplied structure!"
            raise TypeError(err_msg)

        wstru.remove( ligand.key_str )
    
    for lig in ligands:
        lig_cpy = deepcopy(lig)
        wstru.add( lig_cpy, chain_name=lig.key()[0])
        
        relevant_csts:List[StructureConstraint] = list()
        for cst in constraints + stereo_constraints:
            if cst.is_compatible( wstru ):
                cst.change_topology(wstru)
                relevant_csts.append( cst )
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
        opts['ignore_waters'] = False
        opts['include_vrt'] = False
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
        
        dock_ligand(wstru,
                        ligand,
                        relevant_csts,
                        opts,
                        use_qm,
                        cluster_distance,
                        cpu_config
                        )

    stru_oper.update_residues(structure, wstru)
    for cst in constraints + stereo_constraints:
        cst.change_topology( structure )

    sp.save_structure("./snapshots/docked_structure.pdb", structure)
    
    opts.add_script_variable( 'fr_repeats', min_fr_repeats )
    opts.add_script_variable('ramp_constraints', False)
    opts['nstruct'] = 1 #Ooops
    
    mm_minimization(structure, constraints, opts, True, cpu_config)
    sp.save_structure("./snapshots/mm_structure_01.pdb", structure)

    if use_qm:
        qm_minimization_xtb(structure, ligands, constraints, cluster_distance, False, work_dir)
        sp.save_structure("./snapshots/qm_structure_01.pdb", structure)

    #opts.add_script_variable('ramp_constraints', True)
    mm_minimization(structure, constraints, opts, False, cpu_config)
    sp.save_structure("./snapshots/mm_structure_02.pdb", structure)

    if use_qm:
        qm_minimization(structure, ligands, constraints, cluster_distance, False, work_dir, qm_config)
        sp.save_structure("./snapshots/qm_structure_02.pdb", structure)

    translate_structure(structure, start_naming='rosetta')

    if not save_work_dir:
        _LOGGER.info(f"save_work_dir set to False! Deleting {work_dir}")
        fs.safe_rmdir( work_dir )


def mm_minimization(structure:Structure,
                constraints:List[StructureConstraint],
                opts:RosettaOptions,
                bb_flex:bool,
                job_config
                ) -> None:
    """
    """
    
    jump_seles = list()
    for idx, (c1, c2) in enumerate(zip(structure.chains[:-1], structure.chains[1:])):
        jump_seles.append( 
            ('Jump', {'number': str(idx+1), 'setting': str(not(c1.is_polypeptide() and c2.is_polypeptide())) })
        )

    protocol = RosettaScriptsProtocol()
    protocol.add_residue_selector(
		'Index', name="ligand", resnums="%%ligand_idx%%"
    ).add_residue_selector(
		'CloseContact', name="ligand_active_site", residue_selector="ligand", contact_threshold="%%contact_threshold%%"
    ).add_residue_selector(
		'Not', name="not_ligand_active_site", selector="ligand_active_site"
    ).add_scorefunction(
        'ScoreFunction', name='hard_rep', weights='ref2015'        
    ).add_mover(
        'FastRelax', name="frelax", scorefxn="hard_rep", cst_file="%%cst_file%%", 
        repeats="%%fr_repeats%%", 
        ramp_down_constraints="%%ramp_constraints%%", children=[
			('MoveMap', {'name':"full_enzyme", 'bb':"true" if bb_flex else "false", 'chi':"true", 'jump':"false", 'children':[
				('ResidueSelector', {'selector':"ligand_active_site",     'bb': "true" if bb_flex else 'false', 'chi':"true", 'bondangle':"true" }),
				('ResidueSelector', {'selector':"not_ligand_active_site", 'bb': "false", 'chi':"false", 'bondangle':"false"})
                ] + jump_seles })]
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
                    cluster_distance:float=None,
                    job_config=None,
                    ) -> None:
    """TODO(CJ)""" 

    jumps=[]
    for cidx,chain in enumerate(structure.chains):
        if chain.is_polypeptide():
            continue

        if chain.residues[0].key() !=  ligand.key():
            jumps.append(
                ('Jump', {'number': int(cidx), 'setting':'false'})
            )

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
    ).add_residue_selector(
        'ResiduePropertySelector', name='all_ligands', properties='LIGAND'
    ).add_residue_selector(
        'Not', name='not_docked_ligand', selector='ligand'
    ).add_residue_selector(
        'And', name='other_ligands', selectors='not_docked_ligand,all_ligands'
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
        'FastRelax', name="frelax", scorefxn="total_region", cst_file="%%cst_file%%", repeats="%%fr_repeats%%", children=[
            ('MoveMap', {'name':"full_enzyme", 'bb':"true", 'chi':"true", 'jump':"true", 'children':[
                ('ResidueSelector', {'selector':"ligand_active_site", 'bb':"true", 'chi':"true", 'bondangle':"true"}),
				('ResidueSelector', {'selector':"not_ligand_active_site", 'bb':"false", 'chi':"false", 'bondangle':"false"}),
                #('ResidueSelector', {'selector':'other_ligands', 'chi':'false', 'bb':'false', 'bondangle':'false', 'bondlength':'false'}),
                #] + jumps
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
    ligand_key = ligand.key()
    sele_str:str=f"byres all within {cluster_distance*1.5} of (chain {ligand_key[0]} and resi {ligand_key[1]})"
    
    if use_qm:
        protocol.add_scorefunction(
            "ScoreFunction", name='total_region', weights='ref2015'
        ).add_scorefunction(
            "ScoreFunction", name="qm_region", children=[
                ('Reweight', {'scoretype':"orca_qm_energy", 'weight':"1.0" }),
                ('Set', {'orca_path':config['rosetta.ORCA_DIR'] }),
                ('Set', {'orca_processes':"1" }), #TODO(CJ): tunable
                ('Set', {'orca_memory_megabytes':"3000" }), #TODO(CJ): tunable
                ('Set', {'rosetta_orca_bridge_temp_directory':"xtb_temp" }),
                ('Set', {'orca_electron_correlation_treatment':"XTB" }),
                ('Set', {'orca_default_xtb_level': "XTB1"}),
                ('Set', {'clean_rosetta_orca_bridge_temp_directory':"true" }),
                ('Set', {'orca_deduce_charge':"true" })
            ]
        ).add_residue_selector(
            'Index', name='ligand_score', resnums=get_active_site_sele(structure, sele_str, fmt='rosetta', bridge_gaps=True)
        )
        protocol.add_scorefunction(
            'MultiScoreFunction', name='final_scoring', children=[
                ('SimpleCombinationRule',{}),
                ('Region', {'scorefxn':'qm_region', 'residue_selector': 'ligand_score', 'children':[
                    ('CappedBondResolutionRule', {'peptide_nterm_cap':'hydrogen', 'peptide_cterm_cap':'hydrogen'})
                ]}),
                ('Region', {'scorefxn':'total_region', 'children':[
                    ('SimpleBondResolutionRule',)
                ]})
            ]

        )
    else:
        protocol.add_scorefunction(
            'ScoreFunction', name='final_scoring', weights='hard_rep'
        )

    protocol.add_simple_metric(
        'TotalEnergyMetric', name='rosetta_score', scorefxn='final_scoring'
    ).add_protocol(
        metrics='rosetta_score'
    )
    
    ligand_area:float=interface.pymol.get_ligand_area( structure, ligand )
    
    opts.add_script_variable('sasa_cutoff',  int(opts.get_script_variable('sasa_cutoff')*ligand_area))

    opts.add_script_variable('cst_file', 
        interface.rosetta.write_constraint_file( structure, constraints, work_dir=opts['out:path:all'])
    )
    temp_ll = fs.lines_from_file(opts.get_script_variable('cst_file'))
    temp_ct = ''.join(temp_ll).strip()
    assert len(temp_ct), 'no constraints found'
    
    _LOGGER.info("Beginning RosettaLigand geometry sampling step...")
    
    pdb_files = parallel_rs( structure, protocol, opts, job_config, 'nstruct', 5, 'DOCK', './scratch' )

    assert pdb_files

    _parser = PDBParser()
    structures:List[Structure] = list() 
    for pf in pdb_files:
        stru = _parser.get_structure( pf )
        lines = fs.lines_from_file( pf )
        for ll in lines:
            if ll.startswith('rosetta_score'):
                tks = ll.split()
                stru.data['rosetta_score'] = float(tks[-1])
                break
        else:
            assert False, "Missing rosetta_score"
        structures.append( stru )
    
    clusters:List[StructureCluster] = cluster_structures( structures, 'polymer.protein', f"resn {ligand.name}", 1.0 ) #TODO(CJ): update params
    
    ref_stru = sorted(
        clusters,
        key=lambda clust: clust.average_score()
    )[0].lowest_energy_structure()

    stru_oper.update_residues(structure, ref_stru)

    for cst in constraints:
        cst.change_topology(structure)

def get_active_site_sele(structure: Structure, 
                sele_str:str,
                fmt:str='pymol', 
                bridge_gaps:bool=False,
                constraints=None
                ) -> str:
    #TODO(CJ): make this work with constraints
    _LOGGER.info("Analyzing enzyme active site...")
    session = interface.pymol.new_session()
    interface.pymol.load_enzy_htp_stru(session, structure)
    df = interface.pymol.collect(session, 'memory', "chain resi".split(), sele=sele_str)

    result = set()
    for i, row in df.iterrows():
        result.add((row['chain'], int(row['resi'])))

    if constraints is not None:
        for cst in constraints:
            for atom in cst.atoms:
                result.add( atom.parent.key() )

    if bridge_gaps:
        result = sorted(list(result))
        bridge = list()
        for (r1_chain, r1_resi), (r2_chain, r2_resi) in zip(result[:-1], result[1:]):
            if r1_chain != r2_chain:
                continue
            
            r1_resi, r2_resi = int(r1_resi), int(r2_resi)

            if abs(r1_resi - r2_resi) == 2:
                _LOGGER.info(f"Bridging gap between {r1_chain}.{r1_resi} and {r2_chain}.{r2_resi}")
                bridge.append( (r1_chain, r1_resi + 1) )
            
        result.extend( bridge )

    _LOGGER.info(f"Found {len(result)} using selection criteria '{sele_str}'")
    
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
                ligands:List[Ligand],
                constraints:List[StructureConstraint],
                cluster_distance:float,
                freeze_ligands:bool,
                work_dir:str,
                ) -> None:

    sele_str = list()          
    for ll in ligands:
        ligand_key = ll.key()
        sele_str.append(f"(chain {ligand_key[0]} and resi {ligand_key[1]})")
    
    sele_str:str=f"byres all within {cluster_distance} of ( " + " or ".join(sele_str) + " )"
    as_sele:str = get_active_site_sele(
                structure,
                sele_str,
                fmt='pymol',
                bridge_gaps=True,
                constraints=constraints,
                )
    
    constraint_atoms = list()
    for cst in constraints:
        constraint_atoms.extend( cst.atoms )
    
    to_freeze:List[Atom] = list()
    for res in structure.residues:
        if res.is_canonical():
            continue
        
        for atom in res.atoms:
            atom.charge = 0.0

        if res not in ligands:
            for atom in res.atoms:
                if atom not in constraint_atoms and atom.element != 'H':
                    to_freeze.append(atom)                    
        elif freeze_ligands:
            for atom in res.atoms:
                if atom.element != 'H':
                    to_freeze.append(atom)

    translate_structure(structure, start_naming='rosetta')
    es = qm_optimize(structure,
                engine='xtb',
                constraints=constraints + [CartesianFreeze(structure.backbone_atoms() + to_freeze)],
                regions=[as_sele],
                region_methods=[chem.QMLevelOfTheory(basis_set='', method='GFN2', solvent='water', solv_method='ALPB')],
                nterm_cap='H',
                cterm_cap='H',
                parallel_method=None)[0]

    translate_structure(structure, end_naming='rosetta')

def qm_minimization(structure:Structure,
                ligands,
                constraints:List[StructureConstraint],
                cluster_distance:float,
                freeze_ligands:bool,
                work_dir:str,
                job_config
                ) -> None:
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
    sele_str = list()          
    for ll in ligands:
        ligand_key = ll.key()
        sele_str.append(f"(chain {ligand_key[0]} and resi {ligand_key[1]})")
    
    sele_str:str=f"byres all within {cluster_distance} of ( " + " or ".join(sele_str) + " )"

    protocol = RosettaScriptsProtocol()
    options = RosettaOptions()
    protocol.add_residue_selector(
        "Index", name="active_site", resnums=get_active_site_sele(structure, sele_str, 'rosetta', bridge_gaps=False, constraints=constraints)
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
            ( 'Set', {'orca_default_xtb_level': "XTB1"}),
            ( 'Set', {'clean_rosetta_orca_bridge_temp_directory':"false" }),
            ( 'Set', {'orca_deduce_charge':"true" }),
        ]
    ).add_scorefunction(
        'MultiScoreFunction', name="combo_sfxn", children=[
            ('SimplePointChargeCombinationRule', {}),
            ('Region', {'scorefxn':'qm_region', 'residue_selector':'active_site', 'children':[
                ('CappedBondResolutionRule', {'peptide_nterm_cap':'hydrogen', 'peptide_cterm_cap':'hydrogen'})
            ]}),
            ('Region', {'scorefxn':'r15', 'children':[
                ('SimpleBondResolutionRule', {})
            ]}),
        ]
    ).add_mover(
        'OrcaQMGeometryOptimizationMover', name='qm_opt',
            freeze_backbone_atoms='true',
            freeze_CB_atoms='true',
            msfxn_name="combo_sfxn",
            msfxn_freeze_noncommon_atoms="true",
            clean_rosetta_orca_bridge_temp_directory='false',
            geo_opt_max_steps="12", 
            deduce_charge="true",
            orca_memory_megabytes="20000",
            immobilize_h_bond_lengths="true",
            optimization_convergence="NORMALOPT" #TODO(CJ)
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
    options['nstruct'] = 1
    options['ignore_waters'] = False
    options['include_vrt'] = False
    options['keep_input_protonation_state'] = True

    stub_parent: str = os.path.expandvars(
        f"${config['rosetta.ROSETTA3']}/database/chemical/residue_type_sets/fa_standard/residue_types/protonation_states/")
    
    for stub in "GLU_P1 GLU_P2 LYS_D ASP_P1 TYR_D HIS_P ASP_P2".split():
        options.add_extra_res_fa(  f"{stub_parent}/{stub}.params" )

    
    interface.rosetta.parameterize_structure( structure, work_dir )

    for prm in structure.data['rosetta_params']:
        options.add_extra_res_fa( prm )
    
    #TODO(CJ): add in the cartersian freezing here
    constraint_atoms = list()
    for cst in constraints:
        constraint_atoms.extend( cst.atoms )

    to_freeze:List[Atom] = list()
    for res in structure.residues:
        if res.is_canonical(): #TODO(CJ): fix for ASP101
            continue
        
        for atom in res.atoms:
            atom.charge = 0.0

        if res not in ligands:
            for atom in res.atoms:
                if atom not in constraint_atoms and atom.element != 'H':
                    to_freeze.append(atom)                    

        elif freeze_ligands:
            for atom in res.atoms:
                if atom.element != 'H':
                    to_freeze.append(atom)


    print(to_freeze)
    cf = []
    if to_freeze:
        cf = [CartesianFreeze(to_freeze)]
    options.add_script_variable('cst_file', interface.rosetta.write_constraint_file(structure, constraints + cf, 'ORCA_FROZEN', work_dir=work_dir))
    
    options['run:jran'] = 1
    options['out:path:all'] = './'

    infile = parallel_rs( structure, protocol, options, job_config, 'nstruct', 1, 'RQM', './scratch' )[0]
    _parser = PDBParser()
    
    ref_stru = _parser.get_structure(infile)
    stru_oper.update_residues(structure, ref_stru)
    
    for cst in constraints:
        cst.change_topology( structure )


def parallel_rs( 
    structure, 
    protocol, 
    opts, 
    job_config, 
    arr_var, 
    arr_size, 
    arr_prefix, 
    work_dir=None ):

    if work_dir is None:
        work_dir = './'

    orig_rng = opts['run:jran']
    arr_orig = opts[arr_var]
    orig_out_path = opts['out:path:all']
    arrays = [arr_size]*int(arr_orig/arr_size) 
    
    if sum(arrays) != arr_orig:
        arrays.append(  arr_orig%arr_size )
    
    assert sum(arrays) == arr_orig 

    jobs, eggs = list(), list()
    for aidx, arr in enumerate(arrays):
        opts['run:jran'] = aidx
        opts[arr_var] = arr
        opts['out:path:all'] =  str(Path(f"{work_dir}/{arr_prefix}_{aidx:02d}/").absolute())
        fs.safe_rmdir( opts['out:path:all'] )
        engine = RosettaScriptsEngine(
                interface.rosetta, 
                protocol, 
                opts, job_config, opts['out:path:all'])
        (job, egg) = engine.make_job(structure)
        jobs.append( [job] )
        eggs.append( egg )

    jm.ClusterJob.wait_to_2d_array_end(jobs, 60)

    pdb_files = list()
    for ee in eggs:
        pdb_files.extend( ee.pdb_files )
    
    opts[arr_var] = arr_orig 
    opts['out:path:all'] = orig_out_path
    opts['run:jran'] = orig_rng

    result = list()
    for pf in pdb_files:
        if Path(pf).exists():
            result.append( pf )

    return result        

