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
                   ligands:List[Ligand],
                   constraints:List[StructureConstraint]=None,
                   stereo_constraints: List[StructureConstraint]=None,
                   use_qm:bool=True,
                   dock_opts:Dict=None,
                   qm_sele:str=None,
                   qm_freeze_sele:str=None,
                   work_dir:str=None,
                   save_work_dir: bool=True,
                   save_snapshots:bool=True,
                   cpu_config = None,
                   qm_config = None
                   ) -> None:
    """
    """
    if work_dir is None:
        work_dir = config["system.SCRATCH_DIR"]

    if stereo_constraints is None:
        stereo_constraints = []

    fs.safe_mkdir(work_dir)
    if save_snapshots:
        fs.safe_mkdir("./snapshots")
    
    sp = PDBParser()

    if save_snapshots:
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

        
        opts=create_rosetta_opts(structure, ligand, dock_opts)
        
        dock_ligand(wstru, ligand, relevant_csts, opts, use_qm, qm_sele, cpu_config)

    stru_oper.update_residues(structure, wstru)
    for cst in constraints + stereo_constraints:
        cst.change_topology( structure )

    if save_snapshots:
        sp.save_structure("./snapshots/docked_structure.pdb", structure)
    
    mm_minimization(structure, ligands, constraints, qm_sele, opts, cpu_config)

    if save_snapshots:
        sp.save_structure("./snapshots/mm_structure_01.pdb", structure)

    if use_qm:
        qm_minimization(structure, ligands, constraints, qm_sele, qm_freeze_sele, work_dir, qm_config)
        if save_snapshots:
            sp.save_structure("./snapshots/qm_structure_01.pdb", structure)

    translate_structure(structure, start_naming='rosetta')

    if not save_work_dir:
        _LOGGER.info(f"save_work_dir set to False! Deleting {work_dir}")
        fs.safe_rmdir( work_dir )


def mm_minimization(structure:Structure,
                ligands:Ligand,
                constraints:List[StructureConstraint],
                qm_sele:str,
                opts:RosettaOptions,
                job_config
                ) -> None:
    """
    """
  
    opts.add_script_variable( 'fr_repeats', opts.get_script_variable('min_fr_repeats') )
    opts['nstruct'] = 1 #Ooops

    ligand_keys=[r.key() for r in ligands]

    jump_seles = list()
    for idx, ch in enumerate(structure.chains):
        jump_seles.append( 
            ('Jump', {
                'number': str(idx), 
                'setting': (not ch.is_polypeptide() ) and (ch.residues[0].key() in ligand_keys)
        }))
    
    qm_protein_resnums:str=get_active_site_sele(structure, f"({qm_sele}) and polymer.protein", fmt='rosetta')

    protocol = RosettaScriptsProtocol()
    protocol.add_residue_selector(
		'Index', name="asite", resnums=qm_protein_resnums
    ).add_scorefunction(
        'ScoreFunction', name='hard_rep', weights='ligand', children=[
            ('Reweight', {'scoretype':"fa_intra_rep", 'weight':"0.004"}),
            ('Reweight', {'scoretype':"fa_elec", 'weight':"0.42"}),
            ('Reweight', {'scoretype':"hbond_bb_sc", 'weight':"1.3"}),
            ('Reweight', {'scoretype':"hbond_sc", 'weight':"1.3"}),
            ('Reweight', {'scoretype':"rama", 'weight':"0.2"})
        ]
    ).add_mover(
        'FastRelax', name="frelax", scorefxn="hard_rep", cst_file="%%cst_file%%", 
        repeats="%%fr_repeats%%", children=[
			('MoveMap', {'name':"full_enzyme", 'bb':"false", 'chi':"false", 'jump':"false", 'children':[
				('ResidueSelector', {'selector':"asite", 'bb': "true" , 'chi':"true", 'bondangle':"false" }),
                ] + jump_seles })]
    ).add_protocol(
        mover_name='frelax'
    )
    opts.add_script_variable('cst_file', 
        interface.rosetta.write_constraint_file( structure, constraints, work_dir=opts['out:path:all'])
    )
    
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
                    qm_sele:str,
                    job_config=None,
                    ) -> None:
    """TODO(CJ)"""

    jumps=list()
    for cidx,chain in enumerate(structure.chains):
        if chain.is_polypeptide():
            continue

        if chain.residues[0].key() == ligand.key():
            jumps.append(('Jump', {'number': int(cidx), 'setting':'true'}))

    opts['qsar:grid_dir'] = str( Path(opts['out:path:all'] + "/qsar_grids/").absolute())
    fs.safe_mkdir(opts['qsar:grid_dir'])

    qm_resnums:str=get_active_site_sele(structure, qm_sele, fmt='rosetta')
    qm_protein_resnums:str=get_active_site_sele(structure, f"({qm_sele}) and polymer.protein", fmt='rosetta')

    protocol = RosettaScriptsProtocol()
    protocol.add_scoring_grid(
        ligand_chain='%%ligand_chain%%', width='%%grid_width%%', name='grid',children=[
            ('ClassicGrid', {'grid_name':'classic', 'weight':'1.0'}),
            ('HbaGrid', {'grid_name':"hba", 'weight':"2.0"}),
            ('HbdGrid', {'grid_name':"hbd", 'weight':"2.0"})]
    ).add_residue_selector(
        'Index', name='ligand', resnums='%%ligand_idx%%'
    ).add_residue_selector(
        'Index', name="asite", resnums=qm_resnums
    ).add_residue_selector(
        'Index', name="asite_protein", resnums=qm_protein_resnums
    ).add_scorefunction(
        'ScoreFunction', name='ligand_soft_rep', weights='ligand_soft_rep', children=[
            ('Reweight', {'scoretype':"fa_elec", 'weight':"0.42"}),
            ('Reweight', {'scoretype':"hbond_bb_sc", 'weight':"1.3"}),
            ('Reweight', {'scoretype':"hbond_sc", 'weight':"1.3"}),
            ('Reweight', {'scoretype':"rama", 'weight':"0.2"})
        ]
    ).add_scorefunction(
        'ScoreFunction', name='hard_rep', weights='ligand', children=[
            ('Reweight', {'scoretype':"fa_intra_rep", 'weight':"0.004"}),
            ('Reweight', {'scoretype':"fa_elec", 'weight':"0.42"}),
            ('Reweight', {'scoretype':"hbond_bb_sc", 'weight':"1.3"}),
            ('Reweight', {'scoretype':"hbond_sc", 'weight':"1.3"}),
            ('Reweight', {'scoretype':"rama", 'weight':"0.2"})
        ]
    ).add_scorefunction(
        'ScoreFunction', name='total_region'
    ).add_constraint_generator(
        'FileConstraintGenerator', name='add_cst', filename='%%cst_file%%'
    ).add_mover(
        'Transform', name='dock', chain='%%ligand_chain%%', box_size='%%box_size%%', 
            move_distance='%%move_distance%%',  angle="%%transform_angle%%", cycles="%%transform_cycles%%",
            repeats="%%transform_repeats%%", temperature="%%transform_temperature%%", grid_set="grid", 
            use_constraints="true", cst_fa_file="%%cst_file%%"
    ).add_mover(
        'FastRelax', name="frelax", scorefxn="hard_rep", cst_file="%%cst_file%%", repeats="%%fr_repeats%%", children=[
            ('MoveMap', {'name':"full_enzyme", 'bb':"false", 'chi':"false", 'jump':"false", 'children':[
                ('ResidueSelector', {'selector':'asite_protein', 'bb':'false', 'chi':'true', }),
                ] + jumps}
            )]
    ).add_simple_metric(
        'PerResidueClashMetric', name='clash', residue_selector='ligand', residue_selector2='asite'
    ).add_simple_metric(
        'SasaMetric', name="sasa", residue_selector='ligand'
    ).add_simple_metric(
        'TotalEnergyMetric', name='rosetta_score', scorefxn='final_scoring'
    ).add_filter(
        'ConstraintScore', name="cst_filter", constraint_generators="add_cst", threshold="1000000" 
    ).add_protocol(
        mover_name="dock"
    ).add_protocol(
        mover_name="frelax"    
    ).add_protocol(
        filter="cst_filter"
    ).add_protocol(
        metrics="rosetta_score,clash,sasa"                           
    )
    
    if use_qm:
        protocol.add_scorefunction(
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
        ).add_scorefunction(
            'MultiScoreFunction', name='final_scoring', children=[
                ('SimpleCombinationRule',{}),
                ('Region', {'scorefxn':'qm_region', 'residue_selector': 'asite', 'children':[
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

    
    ligand_area:float=interface.pymol.get_ligand_area( structure, ligand )
    
    opts.add_script_variable('sasa_cutoff',  int(opts.get_script_variable('sasa_cutoff')*ligand_area))

    opts.add_script_variable('cst_file', 
        interface.rosetta.write_constraint_file( structure, constraints, work_dir=opts['out:path:all'])
    )
    temp_ll = fs.lines_from_file(opts.get_script_variable('cst_file'))
    temp_ct = ''.join(temp_ll).strip()
    assert len(temp_ct), 'no constraints found'
    
    _LOGGER.info("Beginning RosettaLigand geometry sampling step...")
    
    pdb_files = parallel_rs( structure, protocol, opts, job_config, 'nstruct', 10, 'DOCK', './scratch' )

    assert pdb_files

    _parser = PDBParser()
    structures:List[Structure] = list() 
    for pf in pdb_files:
        variables=collect_data(pf, 'clash cst_filter rosetta_score sasa'.split())
        if variables['clash'] > opts.get_script_variable('clash_cutoff'):
            continue
        if variables['cst_filter'] > opts.get_script_variable('cst_cutoff'):
            continue
        if variables['sasa'] > opts.get_script_variable('sasa_cutoff'):
            continue
        stru = _parser.get_structure( pf )
        stru.data['rosetta_score'] = variables['rosetta_score']
        structures.append( stru )
    
    clusters:List[StructureCluster] = cluster_structures( structures, 'polymer.protein', f"resn {ligand.name}", 1.0 ) #TODO(CJ): update params
    
    ref_stru = sorted(
        clusters,
        key=lambda clust: clust.average_score()
    )[0].lowest_energy_structure()

    stru_oper.update_residues(structure, ref_stru)

    for cst in constraints:
        cst.change_topology(structure)

def collect_data(fname, variables):

    result={}
    for vv in variables:
        result[vv]=None
    lines = fs.lines_from_file( fname )
    for ll in lines:
        for kk in result.keys():
            if not ll.startswith(kk):
                continue
            tks=ll.split()
            if kk==tks[0] or (kk=='clash' and tks[0].find('clash_') != -1):
                result[kk]=float(tks[1])

    return result
    

def get_active_site_sele(structure: Structure, 
                sele_str:str,
                fmt:str='pymol', 
                ) -> str:
    #TODO(CJ): make this work with constraints
    _LOGGER.info("Analyzing enzyme active site...")
    session = interface.pymol.new_session()
    interface.pymol.load_enzy_htp_stru(session, structure)
    df = interface.pymol.collect(session, 'memory', "chain resi".split(), sele=sele_str)

    result = set()
    for i, row in df.iterrows():
        result.add((row['chain'], int(row['resi'])))

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

def qm_minimization(structure:Structure,
                ligands:List[Ligand],
                constraints:List[StructureConstraint],
                qm_sele:str,
                qm_freeze_sele:str,
                work_dir:str,
                cpu_config:Dict
                ) -> None:

    asite_sele=get_active_site_sele(structure, qm_sele)
    session=interface.pymol.new_session()
    (sele, session)=interface.pymol.load_enzy_htp_stru(session, structure)
    df=interface.pymol.collect(
        session, 
        'memory', 
        'chain resi name'.split(), 
        sele=f"(not elem H) and ({qm_freeze_sele}) and ({asite_sele})"
    )
    to_freeze=list()
    for i, row in df.iterrows():
        to_freeze.append(
            structure.get(f'{row.chain}.{row.resi}.{row["name"]}')
        )
    config['xtb.N_PROC'] = 8

    translate_structure(structure, start_naming='rosetta')
    es = qm_optimize(structure,
                engine='xtb',
                constraints=constraints + [CartesianFreeze(structure.backbone_atoms() + to_freeze)],
                regions=[qm_sele],
                region_methods=[chem.QMLevelOfTheory(basis_set='', method='GFN1', solvent='water', solv_method='ALPB')],
                nterm_cap='H',
                cterm_cap='H',
                cluster_job_config=cpu_config,
                job_check_period=20
                )[0]
    translate_structure(structure, end_naming='rosetta')


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
        opts['run:jran'] = aidx+orig_rng
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


def create_rosetta_opts(
                    structure,
                    ligand,
                    dock_opts):

    if not dock_opts:
        dock_opts = {}

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
    opts['run:jran'] = dock_opts.get('rng_seed', 1996)
    opts['nstruct'] = dock_opts.get('n_struct', 100)
    opts['out:path:all'] = dock_opts.get('work_dir', config['system.SCRATCH_DIR'])
    
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
    opts.add_script_variable( 'grid_width', dock_opts.get('grid_width', 50.0) )
    opts.add_script_variable( 'clash_cutoff', dock_opts.get('clash_cutoff', 3) )
    opts.add_script_variable( 'sasa_cutoff', dock_opts.get('max_sasa_ratio', 0.60) ) 
    opts.add_script_variable( 'cst_cutoff', dock_opts.get('cst_energy', 1000) )
    opts.add_script_variable( 'box_size', dock_opts.get('box_size', 50.0) )
    opts.add_script_variable( 'move_distance', dock_opts.get('move_distance', 50.0) )
    opts.add_script_variable( 'transform_angle', dock_opts.get('transform_angle', 360) )
    opts.add_script_variable( 'transform_cycles', dock_opts.get('transform_cycles', 1000) )
    opts.add_script_variable( 'transform_repeats', dock_opts.get('transform_repeats', 3) )
    opts.add_script_variable( 'transform_temperature', dock_opts.get('transform_temperature', 5) )
    opts.add_script_variable( 'fr_repeats', dock_opts.get('docking_fr_repeats', 1) )
    opts.add_script_variable( 'min_fr_repeats', dock_opts.get('min_fr_repeats', 3) )
    
    return opts
