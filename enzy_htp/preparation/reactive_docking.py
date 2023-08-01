"""TODO(CJ)

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-07-28
"""
import os
from pathlib import Path
from typing import List
from copy import deepcopy

import pandas as pd

from enzy_htp import interface

from enzy_htp._interface.rosetta_interface import RosettaCst


from enzy_htp import config
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
                    n_struct:int=None
                    ) -> List[str]:
    """

    Args:

    Returns:

    """

    if not work_dir:
        work_dir = config["system.SCRATCH_DIR"] 

    fs.safe_mkdir( work_dir )

    param_files:List[str] = _prepare_reactants(reactants, reactant_conformers)

    docked_complexes:List[str] = list()
   
    if mutations and constraints:
    #TODO(CJ): make this a prettier function
        kept = []        
        for mut in mutations:
            constrained_mutation=False
            for cst in constraints:
                if mutation_cst_overlap(mut,cst):
                    pass
                    constrained_mutation=True
                    #TODO(CJ): too tired to do this rn
                    exit( 0 )
            if not constrained_mutation:
                kept.append( mut )
        mutations = kept


    (start_pdb, cst_file) = _integrate_csts( structure, constraints, work_dir )

    xml_file:str =  _create_xml(f"{work_dir}/_script.xml", "Y Z".split()) #TODO(CJ): need to get the chain names

    options_file:str = make_options_file(start_pdb, xml_file, param_files, cst_file, work_dir) 
    
    df:pd.DataFrame = _docking_run(options_file)
    
    _evaluate_csts(df, constraints)

    #TODO(CJ): this is where we do the docking run

    #TODO(CJ): here is where we do QM ranking

    

    if mutations:
        for mut in mutations:
            print(mut)
            mutated = mm.mutate_stru( structure, [mut], engine="rosetta" )
            outfile:str = f"{work_dir}/{mut}.pdb"
            print(outfile)
            #TODO(CJ): need to save the file and add to the overall 
                        

    return docked_complexes


def mutation_cst_overlap(mut, cst) -> bool:
    #TODO(CJ): check if the constraint contains a mutation but if it gets mutated back to itself 
    return cst.contains(mut.chain_id, mut.res_idx)


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
            {'parent':'MOVERS', 'tag':'Transform', 'name':f'transform_{cc.lower()}','chain':f'{cc.upper()}',  'box_size':'5', 'move_distance':'0.1', 'angle':'20', 'cycles':'1000', 'repeats':'1', 'temperature':'5'},
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


def make_options_file(pdb_file:str, xml_file:str, param_files:List[str], cst_file:str, work_dir:str) -> str:
    """ """
    content: List[str] = [
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
        f"   -nstruct 5", 
        "   -overwrite",
        "   -path",
        f"       -all './complexes'",
    ])
    
    fname = Path(work_dir) / "options.txt"
    
    fs.write_lines(fname, content)
    
    fs.safe_rmdir(f"{work_dir}/complexes/") 
    fs.safe_mkdir(f"{work_dir}/complexes/")
    
    score_file: str = f"{work_dir}/complexes/score.sc"
    option_file = fname.absolute()
    
    fs.safe_rm(score_file)
    
    
    return str( fname )        
    
    
    
def _docking_run(options:str) -> pd.DataFrame:
    """ """
    
    opt_path = Path(options)
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


def _evaluate_csts(df:pd.DataFrame, csts:List[RosettaCst]) -> pd.DataFrame:
    """ """
    cst_diff = [] 
    for i, row in df.iterrows():
        total = 0.0
        for cst in csts:
            for tl in cst.evaluate( row.description ):
                if tl > 1.0:
                    total += tl
        cst_diff.append( total )

    df['cst_diff'] = cst_diff


def _prepare_reactants(reactants, reactant_conformers) -> List[str]:
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
    

