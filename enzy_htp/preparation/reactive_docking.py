"""Driver for the reactive docking functionality available in enzy_htp. The only function that that should be 
called is dock_reactants(). All others are implementation functions that SHOULD NOT be used. 

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-07-28
"""
import os
from pathlib import Path
from collections import defaultdict
from typing import List, Tuple, Set
from copy import deepcopy

import numpy as np
import pandas as pd

from enzy_htp import interface, config, _LOGGER, binding_energy
from enzy_htp._interface.rosetta_interface import RosettaCst


import enzy_htp.chemical as chem
from enzy_htp import mutation as mm 
from enzy_htp.structure import PDBParser, Structure
from enzy_htp.core import file_system as fs

def dock_reactants( structure:Structure, 
                    reactants:List[str], 
                    constraints:List[RosettaCst] = None, 
                    mutations:List[List[mm.Mutation]] = None,
                    work_dir:str=None,
                    save_work_dir:bool=True,
                    reactant_conformers:List[str]=None,
                    rng_seed:int=1996,
                    n_struct:int=1000,
                    
                    cst_cutoff:int=10,
                    clash_distance:float=2.0,
                    clash_cutoff:int=1,
                    binding_cutoff:int=20,
                    cluster_distance:float=2.5,
                    
                    use_cache:bool=True
                    ) -> List[str]:
    """TODO"""

    if not work_dir:
        work_dir = config["system.SCRATCH_DIR"] 
    
    fs.safe_mkdir( work_dir )

    _place_reactants( structure, reactants )

    (param_files, charge_mapper) = _prepare_reactants(reactants, reactant_conformers)

    docked_complexes:List[str] = list()
   
    resolved = _resolve_mut_cst_overlap(mutations, constraints)
    
    if resolved:
        ((resolved_mutations,resolved_csts),mutants)=resolved
        resolved_complexes = dock_reactants( structure, reactants, resolved_csts, resolved_mutations,
                work_dir,
                save_work_dir,
                reactant_conformers,
                rng_seed,
                n_struct,
                cst_cutoff,
                clash_distance,
                clash_cutoff,
                binding_cutoff,
                cluster_distance,
                use_cache
                ) 

        for rc in resolved_complexes:
            if rc.find('WT') != -1:
                continue
                docked_complexes.append( rc )

    (start_pdb, cst_file) = _integrate_csts( structure, constraints, work_dir, use_cache )

    chain_names:List[str] = _ligand_chain_names( start_pdb )

    xml_file:str =  _create_xml(work_dir, chain_names, use_cache)

    options_file:str = make_options_file(start_pdb, xml_file, param_files, work_dir, rng_seed, n_struct, use_cache, cst_file) 
    
    df:pd.DataFrame = _docking_run(options_file, use_cache)
    
    _evaluate_csts(df, constraints, cst_cutoff)

    _evaluate_clashes(df, clash_distance, clash_cutoff)

    _evaluate_binding(df, start_pdb, binding_cutoff)

    _evaluate_qm(df, start_pdb, charge_mapper, cluster_distance) 

    (selected, wt_fpath) = _select_complex(df, work_dir)

    docked_complexes.append( wt_fpath ) 

    e_flags = []
    for pf in param_files:
        e_flags.extend(['-extra_res_fa', str(Path(pf).absolute())])

    if mutations:
        docked_complexes.extend(
            _mutate_result(selected, mutations, work_dir, e_flags )
        )

    return docked_complexes


def _resolve_mut_cst_overlap(mutations:List[List[mm.Mutation]], constraints:List[RosettaCst]) -> bool:
    """Checks if a mutation and constraint overlap. Note: will return False if a constrained
    residue is mutated back to itself.

    Args:
        mutations: A list(list()) of Mutation() objects to check.
        constraints: A list() of RosettaCst() objects to check.

    Returns:
        Whether the system will try to mutate a constrained residue.
    """
    if not mutations or not constraints:
        return False

    resolved_inputs = list()
    resolved_mutations = list()
    
    for mut_set in mutations:
        
        cst_cpy = []
        overlap = False
        for mut in mut_set:
            for cst in constraints:
                if cst.contains(mut.chain_id, mut.res_idx) and mut.orig != mut.target:
                    overlap = True
                    pass
                else:
                    cst_cpy.append( deepcopy( cst ) )
        
        if overlap:
            resolved_inputs.append((
                [deepcopy(ms) for ms in mut_set],
                cst_cpy
            ))
        else:
            resolved_mutations.append(
                [deepcopy(ms) for ms in mut_set]
            )

    
    if not resolved_inputs:
        return False
    else:
        return (resolved_inputs, resolved_mutations)


def _parameterize_reactant(reactant:str, reactant_name:str, conformers:str=None) -> Tuple[str,int]:
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
        exit( 1 )


    session = interface.pymol.new_session()
    (param_file, pdb_file) = interface.rosetta.parameterize_ligand(reactant,reactant_name, conformers=conformers)
    fs.safe_rm(pdb_file)
    
    sdf_file:str = interface.pymol.convert(session, reactant, new_ext=".sdf")
    fcharge:int = interface.bcl.calculate_formal_charge( sdf_file )

    fs.safe_rm( sdf_file )

    lines:List[str] = fs.lines_from_file(param_file)

    lines.append(
        f"NET_FORMAL_CHARGE {int(fcharge)}"
    )

    fs.safe_rm( param_file )
    fs.write_lines( param_file, lines )

    return (param_file, fcharge)
       

def _integrate_csts( stru:Structure, csts:List[RosettaCst], work_dir:str, use_cache:bool ) -> Tuple[str,str]:
    """
    
    Args:
        stru:
        csts:
        work_dir:

    Returns:
        A Tuple() with the format (pdb file, .cst file).

    """
    parser = PDBParser()
    file_str = parser.get_file_str(stru, if_renumber=False, if_fix_atomname=False)
    
    pdb_content:List[str] = ["HEADER                                            xx-MMM-xx"]
    cst_content:List[str] = list()
    for cidx,cst in enumerate(csts):
        pdb_content.append( cst.create_pdb_line( cidx+1 ) )
        cst_content.extend( cst.create_cst_lines() ) 

    pdb_file:str = f"{work_dir}/start.pdb"
    cst_file:str = f"{work_dir}/rdock.cst"

    if not use_cache:
        fs.safe_rm( pdb_file )
        fs.safe_rm( cst_file )

    if not Path(pdb_file).exists():
        fs.write_lines( pdb_file, pdb_content + file_str.splitlines())

    if not Path(cst_file).exists():
        fs.write_lines( cst_file, cst_content)

    return (pdb_file, cst_file)


def _create_xml(work_dir:str, chains:List[str], use_cache:bool) -> str:
    """ """
    if not use_cache:
        fs.safe_rm(fname)

    fname:str=f"{work_dir}/__script.xml"

    if Path(fname).exists() and use_cache:
        return fname
    #TODO(CJ): do some kind of validation against the chain name
    elements:List[Dict] = [
        {'parent':'SCOREFXNS', 'tag':'ScoreFunction', 'name':'ligand_soft_rep', 'weights':'ligand_soft_rep'},
        {'parent':'SCOREFXNS', 'tag':'ScoreFunction', 'name':'hard_rep', 'weights':'ligand'},
        {'parent':'SCOREFXNS.ScoreFunction', 'tag':'Reweight', 'scoretype':'coordinate_constraint', 'weight':'1.0'},
        {'parent':'SCOREFXNS.ScoreFunction', 'tag':'Reweight', 'scoretype':'atom_pair_constraint', 'weight':'1.0'},
        {'parent':'SCOREFXNS.ScoreFunction', 'tag':'Reweight', 'scoretype':'angle_constraint', 'weight':'1.0'},
        {'parent':'SCOREFXNS.ScoreFunction', 'tag':'Reweight', 'scoretype':'dihedral_constraint', 'weight':'1.0'},
        {'parent':'SCOREFXNS.ScoreFunction', 'tag':'Reweight', 'scoretype':'chainbreak', 'weight':'1.0'},
    ]
    
    for cc in chains:
        elements.extend([
            {'parent':'LIGAND_AREAS', 'tag':'LigandArea', 'name':f'docking_sidechain_{cc.lower()}', 'chain':f'{cc.upper()}', 'cutoff':'6.0', 'add_nbr_radius':'true','all_atom_mode':'false'},
            {'parent':'LIGAND_AREAS', 'tag':'LigandArea', 'name':f'final_sidechain_{cc.lower()}', 'chain':f'{cc.upper()}', 'cutoff':'6.0', 'add_nbr_radius':'true','all_atom_mode':'false'},
            {'parent':'LIGAND_AREAS', 'tag':'LigandArea', 'name':f'final_backbone_{cc.lower()}', 'chain':f'{cc.upper()}', 'cutoff':'7.0', 'add_nbr_radius':'false','all_atom_mode':'true', 'Calpha_restraints':"0.3"},

            {'parent':'INTERFACE_BUILDERS', 'tag':'InterfaceBuilder', 'name':f'side_chain_for_docking_{cc.lower()}', 'ligand_areas':f'docking_sidechain_{cc.lower()}'},
            {'parent':'INTERFACE_BUILDERS', 'tag':'InterfaceBuilder', 'name':f'side_chain_for_final_{cc.lower()}', 'ligand_areas':f'final_sidechain_{cc.lower()}'},
            {'parent':'INTERFACE_BUILDERS', 'tag':'InterfaceBuilder', 'name':f'backbone_{cc.lower()}', 'ligand_areas':f'final_backbone_{cc.lower()}', 'extension_window':'3'},

            {'parent':'MOVEMAP_BUILDERS', 'tag':'MoveMapBuilder', 'name':f'docking_{cc.lower()}', 'sc_interface':f'side_chain_for_docking_{cc.lower()}', 'minimize_water':'false'},
            {'parent':'MOVEMAP_BUILDERS', 'tag':'MoveMapBuilder', 'name':f'final_{cc.lower()}', 'sc_interface':f'side_chain_for_final_{cc.lower()}', 'bb_interface':f'backbone_{cc.lower()}','minimize_water':'false'},


        ])


    elements.extend([
        {'parent':'ROSETTASCRIPTS', 'tag':'SCORINGGRIDS', 'ligand_chain':chains[0].upper(), 'width':'15','append_elements_only':True},
        {'parent':'SCORINGGRIDS', 'tag':'ClassicGrid', 'grid_name':'classic', 'weight':'1.0'},
        {'parent':'MOVERS', 'tag':'AddOrRemoveMatchCsts', 'name':'cstadd', 'cst_instruction':'add_new'}

    ])

    for cc in chains:
        elements.extend([
            {'parent':'MOVERS', 'tag':'Transform', 
                    'name':f'transform_{cc.lower()}','chain':f'{cc.upper()}',  
                    'box_size':'2', 'move_distance':'0.1', 
                    'angle':'45', 'cycles':'500', 'repeats':'5', 'temperature':'25'},
            {'parent':'MOVERS', 'tag':'HighResDocker', 'name':f'high_res_docker_{cc.lower()}', 'cycles':'6', 'repack_every_Nth':'3', 'scorefxn':'ligand_soft_rep', 'movemap_builder':f'docking_{cc.lower()}'},
            {'parent':'MOVERS', 'tag':'FinalMinimizer', 'name':f'final_{cc.lower()}', 'scorefxn':'hard_rep', 'movemap_builder':f'final_{cc.lower()}'}
        ])

    elements.append(
        {'parent':'MOVERS', 'tag':'InterfaceScoreCalculator', 'name':'add_scores', 'chains':','.join(map(lambda ss: ss.upper(), chains)), 'scorefxn':'hard_rep'}
    )


    for cc in chains:
        elements.extend([
            {'parent':'MOVERS', 'tag':'ParsedProtocol', 'name':f'low_res_dock_{cc.lower()}', 'child_nodes':[
                {'tag':'Add', 'mover_name':'cstadd'},
                {'tag':'Add', 'mover_name':f'transform_{cc.lower()}'},
            ]}
        ])

    elements.extend([

        {'parent':'MOVERS', 'tag':'ParsedProtocol', 'name':'high_res_dock', 'child_nodes':
            [deepcopy({'tag':'Add', 'mover_name':f'high_res_docker_{cc.lower()}'}) for cc in chains]
        },
        {'parent':'MOVERS', 'tag':'ParsedProtocol', 'name':'reporting', 'child_nodes':
            [{'tag':'Add', 'mover_name':'add_scores'}]
        },

    ])


    elements.extend(
        [deepcopy({'parent':'PROTOCOLS', 'tag':'Add', 'mover_name':f'low_res_dock_{cc.lower()}'}) for cc in chains] + 
        [{'parent':'PROTOCOLS','tag':'Add', 'mover_name':'high_res_dock'}, {'parent':'PROTOCOLS', 'tag':'Add', 'mover_name':'reporting'}]
    )

    interface.rosetta.write_script(fname, elements)

    return fname


def make_options_file(pdb_file:str,  xml_file:str, param_files:List[str],  work_dir:str, rng_seed:int, n_struct:int, use_cache:bool, cst_file:str=None) -> str:
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
    content: List[str] = [
        "-keep_input_protonation_state",
        "-run:constant_seed",
       f"-run:jran {int(rng_seed)}",
        "-in:file",
       f"    -s '{pdb_file}'",
    ]
    
    for pf in param_files:
        content.append(f"    -extra_res_fa '{Path(pf).absolute()}'")

    stub_parent:str = os.path.expandvars(f"${config['rosetta.ROSETTA3']}/database/chemical/residue_type_sets/fa_standard/residue_types/protonation_states/")
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
        content.extend([
            f"-enzdes:cstfile '{Path(cst_file).absolute()}'"
        ])
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
    
    qsar_grid:str=f"{work_dir}/complexes/qsar_grids/"
    content.append(f"-qsar:grid_dir {qsar_grid}")

    fname = Path(work_dir) / "options.txt"
    score_file: str = f"{work_dir}/complexes/score.sc"
    
    if not use_cache:
        fs.safe_rm(fname)
        fs.safe_rm(score_file)
        fs.safe_rmdir(f"{work_dir}/complexes/") 
        fs.safe_rmdir( qsar_grid )

    fs.safe_mkdir( qsar_grid )
    fs.safe_mkdir(f"{work_dir}/complexes/")

    if not fname.exists():
        fs.write_lines(fname, content)

    option_file = fname.absolute()
    
    return str( fname )        
    
    
    
def _docking_run(option_file:str, use_cache:bool) -> pd.DataFrame:
    """

    Args:

    Returns:

    """

    opt_path = Path(option_file)
    if use_cache: 
        csv_file = Path(str(opt_path.parent / "scores.csv"))
        if csv_file.exists():
            df=pd.read_csv(csv_file)
            df['selected'] = True
            
            for i, row in df.iterrows():
                if not Path(row.description).exists():
                    break
            else:
                return df


    
    start_dir:str = os.getcwd()
    os.chdir(str(opt_path.parent))

    interface.rosetta.run_rosetta_scripts(
        [f"@{opt_path.name}"]
    )

    os.chdir(start_dir)
    
    df:pd.DataFrame=interface.rosetta.parse_score_file(
        str(opt_path.parent / "complexes/score.sc")
    )

    df['description'] = df.apply( lambda row: f"{opt_path.parent}/complexes/{row.description}.pdb", axis=1)
   
    df['selected'] = True

    return df


def _evaluate_csts(df:pd.DataFrame, csts:List[RosettaCst], cst_cutoff:int) -> None:
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
    cst_diff = [] 
    for i, row in df.iterrows():
        total = 0.0
        for cst in csts:
            for tl in cst.evaluate( row.description ):
                if tl > 1.0:
                    total += tl
        cst_diff.append( total )

    df['cst_diff'] = cst_diff
    df['selected'] &= (df['cst_diff'] <= cst_cutoff)


def _prepare_reactants(reactants:List[str], reactant_conformers:List[str]=None) -> List[str]:
    """
    Args:
        reactants:
        reactant_conformers:

    Returns:
        A list() of the .params files for the given reactants.
    """
    #TODO(CJ): error if not enough conformers... then make the conformers
    #TODO(CJ): update function signature
    param_files:List[str] = list()
    charge_mapper:Dict[str,int] = dict()
    for ridx, rct in enumerate(reactants):
        reactant_name:str=f"L{ridx+1:02d}"
        #param_files.append( 
        (pfile, charge) = _parameterize_reactant(
                rct,
                reactant_name,
                conformers=reactant_conformers[ridx] #TODO(CJ): deal with this
            )
        #)
        charge_mapper[reactant_name] = charge
        param_files.append( pfile )

    return (param_files, charge_mapper)

def _dock_mut_cst_overlaps(structure:Structure,
                            reactants:List[str], 
                            constraints:List[RosettaCst], 
                            mutations:List[mm.Mutation],
                            work_dir:str,
                            save_work_dir:bool,
                            reactant_conformers:List[str],
                            rng_seed:int,
                            n_struct:int
                            ) -> List[str]:
    """
    
    Args:

    Returns:

    """
    pass
    #TODO(CJ): implement this


def _mutate_result(structure, mutations:List[List[mm.Mutation]], work_dir :str, extra_flags=None) -> List[str]:
    """ TODO(CJ)"""

    res_names = {}
    for res in structure.residues:
        if res.name.upper() in chem.METAL_CENTER_MAP:
            continue
        if res.name not in chem.THREE_LETTER_AA_MAPPER:
            continue
        res_names[
            f"{res.chain.name}.{res.idx}"
            ] = chem.convert_to_one_letter(res.name)

    if extra_flags is None:
        extra_flags = list()

    result:List[str] = list()
    #TODO(CJ): get this to work for List[List[Mutation]]
    for mut_sets in mutations:
        mutated_outfile=str()
        for mut in mut_sets:
            original = res_names[f"{mut.chain_id}.{mut.res_idx}"]
            mutated_outfile += f"{original}{mut.res_idx}{chem.convert_to_one_letter(mut.target)}_"
        
        mutated_outfile:str = f"{work_dir}/{mutated_outfile[:-1]}.pdb"
        mutant = mm.mutate_stru(structure, mut_sets, engine='rosetta', extra_flags=extra_flags )
        parser = PDBParser()
        file_str = parser.get_file_str(mutant, if_renumber=False, if_fix_atomname=False)
        fs.safe_rm(mutated_outfile)
        fs.write_lines( mutated_outfile, file_str.splitlines() )
        result.append( mutated_outfile )

    return result 


def _ligand_chain_names( start_pdb : str ) -> List[str]:
    """
    Args:

    Returns:
    """
    fs.check_file_exists( start_pdb ) 
    session = interface.pymol.new_session()
    df:pd.DataFrame=interface.pymol.collect(session, start_pdb, "chain resn".split()    )

    result:List[str] = list()

    for i, row in df.iterrows():

        if row.resn in chem.THREE_LETTER_AA_MAPPER:
            continue

        if row.resn.upper() in chem.METAL_CENTER_MAP:
            continue

        result.append(row.chain)


    return list(set(result))


def _select_complex(df : pd.DataFrame, work_dir:str, cutoff:int=20) -> Tuple[Structure,str]:
    """Final selection of the complex that we will use.
    
    Args:
        df:
        work_dir:
        cutoff:

    Returns:
        A Tuple() with the layout (selected structure, selected .pdb file).
    """
    #TODO(CJ): I need to change this so that I am just ranking based on QM energy
    # then I should:
    # 1. save the dataframe  -> still need to do this
    # 2. get the structure
    # 3. get the wt filepath
    cpy = deepcopy( df )
    
    cpy = cpy[cpy.selected].reset_index(drop=True)

    cpy.sort_values(by='qm_energy', inplace=True)

    wt_fpath:str=cpy.description[0] 

    #TODO(CJ): check that the df does  have length of at least 1

    parser = PDBParser()

    result_stru = parser.get_structure( wt_fpath )

    return (result_stru, str(Path(wt_fpath).absolute()) )




def _place_reactants( structure:Structure, reactants:List[str] ) -> bool:
    """TODO(CJ)"""
    
    n_reactants:int=len(reactants)
    reactants_left:Set=set([Path(rr).stem for rr in reactants])
    
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
        exit( 1 )


def _evaluate_clashes(df:pd.DataFrame, clash_distance:float, clash_cutoff:int) -> None:
    """TODO(CJ)"""

    session = interface.pymol.new_session()
    clash_ct:List[int] = list() 
    for i, row in df.iterrows():
        
        if not row.selected:
            clash_ct.append( -1 )
            continue
        
        atoms:pd.DataFrame = interface.pymol.collect(session, row.description, "elem resn vdw x y z".split())
        atoms:pd.DataFrame = atoms[atoms.elem!='H'].reset_index(drop=True)

        ligand_names = list(
            filter( lambda rn: rn[0] == 'L' and rn[1].isnumeric() and rn[2].isnumeric(), atoms.resn.unique())
        )

        other = atoms[~(atoms.resn.isin(ligand_names))].reset_index(drop=True)
        o_xyz = np.transpose(np.array([other.x.to_numpy(), other.y.to_numpy(), other.z.to_numpy()]))
        
        count:int = 0
        for resn in ligand_names:
            ligand = atoms[atoms.resn==resn].reset_index(drop=True)
            l_xyz = np.transpose(np.array([ligand.x.to_numpy(), ligand.y.to_numpy(), ligand.z.to_numpy()]))
    

            for ll in l_xyz:
                if np.min(
                    np.sqrt(np.sum(np.power(o_xyz-ll,2),axis=1))
                    ) <= clash_distance:
                    count += 1
        clash_ct.append( count )
    
    df['clash_ct'] = clash_ct

    df['selected'] = (df.selected)&(df.clash_ct<=clash_cutoff)
    

def _evaluate_binding(df:pd.DataFrame, start_pdb:str, binding_cutoff:int, distance_cutoff:float=4.0) -> None:
    """ """
    session = interface.pymol.new_session()
    res_names:List[str] = interface.pymol.collect(session, start_pdb, "resn".split()).resn.unique()
    
    sele_names:Dict[str, str] = {}
    binding_ddg:Dict[str, List[float]] = {}
    charges:Dict[str,float] = {}


    temp_sdf:str=f"{config['system.SCRATCH_DIR']}/__temp.sdf"
    fs.safe_rm(temp_sdf)


    for rn in res_names:
        if rn.upper() in chem.METAL_CENTER_MAP or rn.upper() in chem.THREE_LETTER_AA_MAPPER:
            continue
        
        atoms:pd.DataFrame = interface.pymol.collect(session, start_pdb,
            "resi chain".split(),
            sele=f"(all within {distance_cutoff:.2f} of resn {rn}) and not resn {rn}"
            )
        
        sele_names[rn] = " or ".join(map(lambda pr: f"(resi {pr[0]} and chain {pr[1]})",set(zip(atoms.resi,atoms.chain))))
        binding_ddg[rn] = list()

        receptor_cluster:str=interface.pymol.create_cluster(session, 
                    start_pdb,
                    sele_names[rn],
                    cap_strategy='CH3',
                    outfile=temp_sdf
                    )
        
        charges[f"{rn}_receptor"] = interface.bcl.calculate_formal_charge( temp_sdf)
        
        probe_cluster:str=interface.pymol.create_cluster(session, 
                    start_pdb,
                    f"resn {rn}",
                    cap_strategy='CH3',
                    outfile=temp_sdf
                    )

        charges[f"{rn}_probe"] = interface.bcl.calculate_formal_charge( temp_sdf )
        
    fs.safe_rm( temp_sdf )
    reactants:List[str] = list(sele_names.keys())

    for i, row in df.iterrows():
        if not row.selected:
            for rn in reactants:
                binding_ddg[rn].append( 1.0 )
            continue

        
        for rn in reactants:
            binding_ddg[rn].append(
                binding_energy(
                    row.description,
                    f"resn {rn}",
                    sele_names[rn],
                    charges[f"{rn}_probe"],
                    charges[f"{rn}_receptor"],
                    work_dir=config["system.SCRATCH_DIR"]

                )
            )
  
    all_negative = np.array([True]*len(df)) 
    binding_log_sum = np.zeros(len(df))     
    for rn, values in binding_ddg.items():
        df[f"{rn}_ddg"] = values
        binding_log_sum += np.exp(df[f"{rn}_ddg"].to_numpy())
        all_negative &= (df[f"{rn}_ddg"] <= 0.0)

    df['all_negative'] = all_negative
    df['binding_log_sum'] = binding_log_sum
    
    if not df.selected.sum():
        return

    bls_cutoff:float=np.percentile(df[df.selected].binding_log_sum, binding_cutoff)
    
    an_cutoff = int((df.selected.sum())*(float(binding_cutoff)/100.0))
    
    if df.all_negative.sum() >= an_cutoff:
        df['selected'] = (df.selected)&(df.all_negative)
    else:
        df['selected'] = (df.selected)&(df.binding_log_sum<=cutoff)


def _evaluate_qm(df:pd.DataFrame, start_pdb:str, charge_mapper, cluster_cutoff:float) -> None:
    """ """
    # steps
    # 1. get system charge
    # 2. run through each cluster and do it
    qm_energy:List[float] = list()
    
    session = interface.pymol.new_session()
    res_keys = interface.pymol.collect(session, start_pdb, "chain resn resi".split())
    res_keys = sorted(set(zip(res_keys.chain, res_keys.resn, res_keys.resi)))

    reactants = set()
    reactant_counter = defaultdict(int)
    
    for (cname, resn, resi) in res_keys:
        if resn.upper() in chem.METAL_CENTER_MAP or resn.upper() in chem.THREE_LETTER_AA_MAPPER:
            continue
        
        reactants.add(resn)

        reactant_counter[resn] += 1


    as_sele:str=" or ".join(map(lambda rn: f"(byres all within {cluster_cutoff} of resn {rn})",reactants))
    
    interface.pymol.general_cmd(session,[('delete','all')])
    res_keys = interface.pymol.collect(session, start_pdb, "chain resn resi".split(), sele=as_sele)
    res_keys = sorted(set(zip(res_keys.chain, res_keys.resn, res_keys.resi)))
    atom_mapper = defaultdict(list)
    atoms = interface.pymol.collect(session, start_pdb, "chain resn resi name".split(), sele=as_sele)
    for i, row in atoms.iterrows():
        atom_mapper[(row.chain, row.resn, row.resi)].append( row['name'] )

    charge = 0
    
    for (name, res_name, resi) in res_keys:
        if res_name in chem.AA_CHARGE_MAPPER:
            charge += chem.AA_CHARGE_MAPPER[res_name]
        elif res_name in charge_mapper:
            charge += charge_mapper[res_name]
        elif res_name.upper() == "MG":
            charge += 2
        else:
            assert False, res_name

        if res_name == 'ASP':
            res_atoms = atom_mapper[(name, res_name, resi)]
            if 'HD2' in res_atoms:
                charge += 1

    for i, row in df.iterrows():
        if not row.selected:
            qm_energy.append( 0  )
            continue


        interface.pymol.general_cmd(session,[('delete','all')])
        interface.pymol.create_cluster(session, row.description, as_sele, outfile='temp.mol', cap_strategy='CH3')
        qm_energy.append( 
            interface.xtb.single_point('temp.mol', charge=charge)
        )


    df['qm_energy'] = qm_energy
