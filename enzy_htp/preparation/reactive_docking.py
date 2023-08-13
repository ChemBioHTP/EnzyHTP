"""TODO(CJ)

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-07-28
"""
import os
from pathlib import Path
from typing import List
from copy import deepcopy

import numpy as np
import pandas as pd

from enzy_htp import interface

from enzy_htp._interface.rosetta_interface import RosettaCst


from enzy_htp import config
import enzy_htp.chemical as chem
from enzy_htp import mutation as mm 
from enzy_htp.structure import PDBParser
from enzy_htp.core import file_system as fs


#TODO(CJ): add an rng seed
def dock_reactants( structure, 
                    reactants:List[str], 
                    constraints:List[RosettaCst] = None, 
                    mutations:List[mm.Mutation] = None,
                    work_dir:str=None,
                    save_work_dir:bool=False,
                    reactant_conformers:str=None,
                    rng_seed:int=1996,
                    n_struct:int=1000,
                    use_cache:bool=True
                    ) -> List[str]:
    """

    Args:
        structure:
        reactants:
        contraints:
        mutations:
        work_dir:
        save_work_dir:
        reactant_conformers:
        rng_seed:
        n_struct:
        use_cache:

    Returns:

    """
    #TODO(CJ): check if the ligands are in the structure
    #TODO(CJ): checks on 
    if not work_dir:
        work_dir = config["system.SCRATCH_DIR"] 

    fs.safe_mkdir( work_dir )

    param_files:List[str] = _prepare_reactants(reactants, reactant_conformers)

    docked_complexes:List[str] = list()
   
    if (mutations and constraints) and mut_cst_overlap(mutations, constraints):
        docked_complexes.extend(
            _dock_mut_cst_overlaps(
                        structure, 
                        reactants, 
                        constraints, 
                        mutations,
                        work_dir,
                        save_work_dir,
                        reactant_conformers,
                        rng_seed,
                        n_struct
                    )
        )
    
    (start_pdb, cst_file) = _integrate_csts( structure, constraints, work_dir )

    chain_names:List[str] = _ligand_chain_names( start_pdb )

    xml_file:str =  _create_xml(f"{work_dir}/_script.xml", chain_names)

    options_file:str = make_options_file(start_pdb, xml_file, param_files, cst_file, work_dir, rng_seed, n_struct) 
    
    df:pd.DataFrame = _docking_run(options_file)
    
    _evaluate_csts(df, constraints)

    df['clash_ct'] = df.apply(lambda row: _clash_count(row.description), axis=1)

    df.to_csv(f"{work_dir}/scores.csv", index=False)
    (selected, wt_fpath) = _select_complex(df, work_dir)

    docked_complexes.append( wt_fpath ) 
    #TODO(CJ): here is where we do QM ranking

    if mutations:
        docked_complexes.extend(
            _mutate_result(selected, mutations, work_dir )
        )

    df.to_csv(f"{work_dir}/scores.csv", index=False)

    return docked_complexes


def mut_cst_overlap(mutations, constraints) -> bool:
    #TODO(CJ): check if the constraint contains a mutation but if it gets mutated back to itself 
    for mut in mutations:
        for cst in constraints:
            if cst.contains(mut.chain_id, mut.res_idx):
                return True
    return False


def _parameterize_reactant(reactant:str, reactant_name:str, conformers:str=None) -> str:
    """ """

    session = interface.pymol.new_session()
    #TODO(CJ): this is where the parameterizing goes
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

    return param_file
       

def _integrate_csts( stru, csts, work_dir ):
    parser = PDBParser()
    file_str = parser.get_file_str(stru, if_renumber=False, if_fix_atomname=False)
    
    pdb_content:List[str] = ["HEADER                                            xx-MMM-xx"]
    cst_content:List[str] = list()
    for cidx,cst in enumerate(csts):
        pdb_content.append( cst.create_pdb_line( cidx+1 ) )
        cst_content.extend( cst.create_cst_lines() ) 

    pdb_file:str = f"{work_dir}/start.pdb"
    cst_file:str = f"{work_dir}/rdock.cst"

    fs.safe_rm( pdb_file )
    fs.safe_rm( cst_file )

    fs.write_lines( pdb_file, pdb_content + file_str.splitlines())
    fs.write_lines( cst_file, cst_content)
    return (pdb_file, cst_file)


def _create_xml(fname:str, chains:List[str]) -> str:
    """ """
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


def make_options_file(pdb_file:str,  xml_file:str, param_files:List[str], cst_file:str, work_dir:str, rng_seed:int, n_struct:int) -> str:
    """ """
    content: List[str] = [
        "-run:constant_seed",
       f"-run:jran {int(rng_seed)}",
        "-in:file",
       f"    -s '{pdb_file}'",
    ]
    
    for pf in param_files:
        content.append(f"    -extra_res_fa '{Path(pf).absolute()}'")

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
        #TODO(CJ): add a warning if there are no .csts
        content.extend([
            f"-enzdes:cstfile '{Path(cst_file).absolute()}'"
        ])
    
    
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
    
    fs.write_lines(fname, content)
    
    fs.safe_rmdir(f"{work_dir}/complexes/") 
    fs.safe_mkdir(f"{work_dir}/complexes/")
   
    fs.safe_rmdir( qsar_grid )
    fs.safe_mkdir( qsar_grid )

    score_file: str = f"{work_dir}/complexes/score.sc"

    option_file = fname.absolute()
    
    fs.safe_rm(score_file)
    
    
    return str( fname )        
    
    
    
def _docking_run(option_file:str) -> pd.DataFrame:
    """

    Args:

    Returns:

    """
    
    opt_path = Path(option_file)
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
    
    return df


def _evaluate_csts(df:pd.DataFrame, csts:List[RosettaCst]) -> None:
    """
    Args:
        df:
        csts:

    Returns:
        Nothing.
    """
    cst_diff = [] 
    for i, row in df.iterrows():
        total = 0.0
        for cst in csts:
            #TODO(CJ): Check if the row.description attribute works
            for tl in cst.evaluate( row.description ):
                if tl > 1.0:
                    total += tl
        cst_diff.append( total )

    df['cst_diff'] = cst_diff


def _prepare_reactants(reactants:List[str], reactant_conformers:List[str]=None) -> List[str]:
    """
    Args:
        reactants:
        reactant_conformers:

    Returns:
        A list() of the .params files for the given reactants.
    """
    #TODO(CJ): error if not enough conformers... then make the conformers
    param_files:List[str] = list()
    for ridx, rct in enumerate(reactants):
        reactant_name:str=f"L{ridx+1:02d}"
        param_files.append( 
            _parameterize_reactant(
                rct,
                reactant_name,
                conformers=reactant_conformers[ridx] #TODO(CJ): deal with this
            )
        )

    return param_files

def _dock_mut_cst_overlaps(structure,
                            reactants:List[str], 
                            constraints:List[RosettaCst], 
                            mutations:List[mm.Mutation],
                            work_dir:str,
                            save_work_dir:bool,
                            reactant_conformers:List[str],
                            rng_seed:int,
                            n_struct:int
                            ) -> List[str]:
    """ """
    pass
    #TODO(CJ): implement this


def _mutate_result(structure, mutations:List[mm.Mutation], work_dir :str) -> List[str]:
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

    result:List[str] = list()
    for mut in mutations:
        original = res_names[f"{mut.chain_id}.{mut.res_idx}"]
        mutated_outfile:str=f"{work_dir}/{original}{mut.res_idx}{chem.convert_to_one_letter(mut.target)}.pdb"
        mutant = mm.mutate_stru(structure, [mut], engine='rosetta')
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

def _clash_count(fname:str, cutoff:float=2.25) -> int:
    """ """
    session = interface.pymol.new_session()
    df = interface.pymol.collect(session, fname, "elem resn vdw x y z".split())
    df:pd.DataFrame = df[df.elem!='H']
    df.reset_index(drop=True,inplace=True)

    ligand_names = list()
    for rn in df.resn.unique():
        if rn[0] == 'L' and rn[1].isnumeric() and rn[2].isnumeric():
            ligand_names.append( rn )

    count = 0
    for resn in ligand_names:
        ligand = df[df.resn==resn].reset_index(drop=True)
        l_xyz = np.transpose(np.array([ligand.x.to_numpy(), ligand.y.to_numpy(), ligand.z.to_numpy()]))
        other = df[df.resn!=resn].reset_index(drop=True)
        o_xyz = np.transpose(np.array([other.x.to_numpy(), other.y.to_numpy(), other.z.to_numpy()]))

        for ll in l_xyz:
            if np.min(np.sqrt(np.sum(np.power(o_xyz-ll,2),axis=1))) < cutoff:
                count += 1

    return count


def _select_complex(df : pd.DataFrame, work_dir:str, cutoff:float=0.20):
    n = max(1, int(cutoff*len(df)))
    df.sort_values(by='total_score', inplace=True)
    e_cutoff = df.total_score.to_list()[n]
    df.sort_values(by='cst_diff', inplace=True)
    cst_cutoff = df.cst_diff.to_list()[n]
    df.sort_values(by='clash_ct', inplace=True)
    clash_cutoff = df.clash_ct.to_list()[n]

    mask = (df.total_score>=e_cutoff)&(df.cst_diff>=cst_cutoff)&(df.clash_ct>=clash_cutoff)
    
    print(sum(mask))

    for i, row in df[mask].iterrows():
        print(row.description)        
    
    #energy_sorted = df.sort_values(by='total_score')
    #energy_sorted = energy_sorted.head(n)
    #cst_sorted = df.sort_values(by='cst_diff')
    #
    exit(0)
    parser = PDBParser()
    selected = parser.get_structure(cst_sorted.description[0])
    outfile = f"{work_dir}/WT.pdb"
    
    file_str = parser.get_file_str(selected, if_renumber=False, if_fix_atomname=False)
    fs.write_lines(outfile, file_str.splitlines())

    return (selected, outfile)

