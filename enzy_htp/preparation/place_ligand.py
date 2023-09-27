"""Provides functionality for placing a ligand into a protein. All functionalty should be accessed through 
the place_ligand() function. No other function should be called directly by users. Allows for placement of 
ligands via one of the below methods:

    + alphafill: Uses the AlphaFill transplant algorithm to identify a location for a given ligand. Works best
    for common bi-molecular co-factors (i.e. SAM, ATP, etc). Assumes PDB style atom naming.
    + mole2: TODO(CJ)

Author: Chris Jurich <chris.jurich@vanderbit.edu>
Date: 2023-09-14
"""

from typing import List, Tuple, Dict

from rdkit import Chem as _rchem
import numpy as np
import pandas as pd
import enzy_htp as eh
from pathlib import Path
from copy import deepcopy
from collections import defaultdict

from enzy_htp import interface, config
from enzy_htp.core import file_system as fs

from enzy_htp.core import _LOGGER

import enzy_htp.chemical as chem
from enzy_htp.structure import PDBParser, Structure
from enzy_htp._interface import RosettaCst, Mole2Cavity

from .align_ligand import align_ligand
from .constrained_docking import constrained_docking

#TODO(CJ): for the RosettaLigand version, set the COM to the COM of the constrained residue
def place_ligand(molfile: str,
                 ligand: str,
                 code: str = None,
                 method: str = 'alphafill',
                 new_res_key: Tuple = None,
                 work_dir: str = None,
                 outfile: str = None,
                 constraints: List[RosettaCst] = None,
                 n_struct:int = 50,
                 use_cache:bool = True,
                 **kwargs) -> Tuple[str, str]:
    """Protocol for placing a ligand into an enzyme,  

    Args:
        molfile:
        ligand:
        code:
        method:
        new_res_key:
        work_dir:
        outfile:
        constraints:
        n_struct: 


    Returns:
        A Tuple[str,str] with the layout (enzyme-ligand complex in .pdb format, exact ligand conformation in .mol2 format).
    """
    #TODO(CJ): check for None values and throw errors if they are given
    if not work_dir:
        work_dir = config['system.SCRATCH_DIR']

    fs.safe_mkdir(work_dir)

    placed_ligand: str = None

    if method == "alphafill":
        placed_ligand = _place_alphafill(molfile, ligand, code=code, work_dir=work_dir, **kwargs)
    elif method == "mole2":
        pass
        placed_ligand = _place_mole2(molfile, ligand, constraints=constraints, new_res_key=new_res_key, work_dir=work_dir, **kwargs)
        #place_ligand = _place_rosetta_ligand(molfile, ligand, new_res_key=new_res_key, work_dir=work_dir, constraints=constraints, n_struct=n_struct, use_cache=use_cache, **kwargs)
    else:
        _LOGGER.error(
            f"The supplied method placement method {method} is not supported. Allowed methods are 'alphafill' and 'mole2'. Exiting..."
        )
        exit(1)

    session = interface.pymol.new_session()

    args = [('load', placed_ligand)]

    if new_res_key is not None:
        args.extend([
            ('alter', 'all', f'chain="{new_res_key[0]}"'),
            ('alter', 'all', f'resi="{new_res_key[1]}"'),
            ('alter', 'all', f'resn="{new_res_key[2]}"'),
            ('alter', 'all', 'segi=""'),
        ])

    if not outfile:
        temp_path = Path(molfile)
        outfile = f"{work_dir}/{temp_path.stem}_placed{temp_path.suffix}"

    args.extend([('save', placed_ligand), ('load', molfile), ('save', outfile)])

    interface.pymol.general_cmd(session, args)

    return (outfile, placed_ligand)

def _place_mole2(molfile, ligand, constraints, new_res_key, work_dir:str, **kwargs) -> str:
    """TODO(CJ)"""
    #TODO(CJ): check how many atoms fit in the molecule
    if not work_dir:
        work_dir = config['system.SCRATCH_DIR']

    relevant_constraints: List[str] = list()
    
    #TODO(CJ): probably need to deepcopy these
    for cc in constraints:
        if cc.contains(new_res_key[0], new_res_key[1]):
            relevant_constraints.append(cc)

    for rc in relevant_constraints:
        rc.remove_constraint('angle_A')
        rc.remove_constraint('angle_B')
        rc.remove_constraint('torsion_A')
        rc.remove_constraint('torsion_B')
        rc.remove_constraint('torsion_AB')

    anchors = list()        
    for rc in relevant_constraints:
        anchors.append(rc.other(new_res_key[0], new_res_key[1]))


    cavities:List[Mole2Cavity] = interface.mole2.identify_cavities(molfile)
    for cc in cavities:
        print(cc, cc.center_of_mass(), cc.volume())

    session = interface.pymol.new_session()
    interface.pymol.general_cmd(session, [('load', molfile)])
    anchor_com = interface.pymol.center_of_mass(session, sele=f"chain {anchors[0][0]} and resi {anchors[0][1]}")
   
    distances = list()
    for cc in cavities:
        distances.append(np.sqrt(np.sum((cc.center_of_mass()-anchor_com)**2)))
   

    target_cat = cavities[np.argmin(np.array(distances))]
    
    target_com = target_cat.center_of_mass()

    interface.pymol.general_cmd(session, [('delete', 'all'), ('load', ligand)])
    reactant_com = interface.pymol.center_of_mass(session)
    shift = target_com - reactant_com

    temp_path = Path(ligand)
    outfile = str(Path(work_dir) / f"{temp_path.stem}_shifted{temp_path.suffix}")
    interface.pymol.general_cmd(session,[
        ('alter_state', '1', '(all)', f"x+={shift[0]}"), 
        ('alter_state', '1', '(all)', f"y+={shift[1]}"), 
        ('alter_state', '1', '(all)', f"z+={shift[2]}"), 
        ('save', outfile)
    ])

    return outfile


def _place_rosetta_ligand(molfile: str, reactant: str, new_res_key: Tuple, work_dir: str, constraints: List[RosettaCst], n_struct:int, use_cache:bool, **kwargs) -> str:
    """Protocol that uses Rosetta's RosettaLigand and enzdes constraint functionality to place a ligand into an enzyme. SHOULD NOT be called
    directly by the user. Instead, call `enzy_htp.preparation.place_ligand.place_ligand()`. Uses the starting file, reactant, new residue key, 
    and constraints to perform a fast docking protocol which does not perform high resolution docking or side-chain packaing. Note that at least
    one of the supplied constraints must contain the ligand to be placed. If it does not, an error will be thrown. 

    Args:
        molfile:
        reactant:
        new_res_key:
        work_dir: 
        constraints:


    Returns:
        
    """
    if not work_dir:
        work_dir = config['system.SCRATCH_DIR']

    relevant_constraints: List[str] = list()

    for cc in constraints:
        if cc.contains(new_res_key[0], new_res_key[1]):
            relevant_constraints.append(cc)

    for rc in relevant_constraints:
        rc.remove_constraint('angle_A')
        rc.remove_constraint('angle_B')
        rc.remove_constraint('torsion_A')
        rc.remove_constraint('torsion_B')
        rc.remove_constraint('torsion_AB')

    print(relevant_constraints)

    #TODO(CJ): find the anchor residue 
    anchors = list()        
    for rc in relevant_constraints:
        anchors.append(rc.other(new_res_key[0], new_res_key[1]))

    shifted_reactant:str = _shift_reactant(molfile, reactant, anchors, work_dir)

    #!move it to near the anchor!!

    if not relevant_constraints:
        _LOGGER.error(
            f"To place a ligand with RosettaLigand, you must supply at least one constraint involving the target ligand. Exiting...")
        exit(1)
    else:
        _LOGGER.info(f"Found {len(relevant_constraints)} constraints that involve the ligand! Continuing...")

    temppath = Path(molfile)
    other = str(temppath.parent / "other.pdb")

    session = interface.pymol.new_session()
    interface.pymol.general_cmd(
        session,
        [
            ('load', shifted_reactant),
            ('alter', 'all', f"chain='{new_res_key[0]}'"),
            ('alter', 'all', f"resi='{new_res_key[1]}'"),
            ('alter', 'all', f"resn='{new_res_key[2]}'"),
            ('alter', 'all', f"segi=''"),
            #('load', molfile),
            ('save', shifted_reactant),
            ('delete', 'all')
        ])
    interface.pymol.general_cmd(session, [('load', molfile)])

    df: pd.DataFrame = interface.pymol.collect(session, 'memory', 'chain resn resi'.split())
    need_params: List[str] = list()

    for i, row in df.iterrows():
        key = (row.chain, row.resn, row.resi)

        if key in need_params:
            continue

        if key[1].upper() in chem.METAL_CENTER_MAP:
            continue

        if key[1].upper() in chem.THREE_LETTER_AA_MAPPER:
            continue

        need_params.append(key)

    param_files: List[str] = list()
    for needs in need_params:
        temp_file: str = f"{work_dir}/temp.mol2"
        interface.pymol.general_cmd(session, [("save", temp_file, f"chain {needs[0]} and resn {needs[1]} and resi {needs[2]}")])
        param_files.append(interface.rosetta.parameterize_ligand(temp_file, needs[1])[0])

    param_files.append(interface.rosetta.parameterize_ligand(reactant, new_res_key[2])[0])
    parser = eh.PDBParser()

    interface.pymol.general_cmd(session, [('load', shifted_reactant), ('save', other)])

    structure = parser.get_structure(other)
    #TODO(CJ): make this less error prone
    #TODO(CJ): need to look through and find all the parts where there are
    df:pd.DataFrame=constrained_docking(structure, [reactant], param_files, relevant_constraints, work_dir=work_dir, minimize=False, n_struct=n_struct, use_cache=use_cache)
    cst_diff:List[float] = list()
    for i, row in df.iterrows():
        total:float = 0.0
        for rc in relevant_constraints:
            total += sum(rc.evaluate(row.description))
        cst_diff.append( total )
    df['cst_diff'] = cst_diff
    print(cst_diff)
    print(df)
    #what needs to happen here?
    # 1. figure out which constraints involve the reactant we care about
    # 2. create the input files that we'll need
    # 3. run the rosetaLigand run
    # 4. select placement with best constraint score
    assert False

def _shift_reactant(molfile:str, reactant:str, anchors:List[str], work_dir:str) -> str:
    """
    Args:

    Returns:
        The shifted
    """
    #TODO(CJ): some code about selecting the right anchor constraint        
    session = interface.pymol.new_session()
    interface.pymol.general_cmd(session, [('load', molfile)])
    anchor_com = interface.pymol.center_of_mass(session, sele=f"chain {anchors[0][0]} and resi {anchors[0][1]}")
    #TODO(CJ): add error if the com is zero
    interface.pymol.general_cmd(session, [('delete', 'all'), ('load', reactant)])
    reactant_com = interface.pymol.center_of_mass(session)
    shift = anchor_com - reactant_com

    temp_path = Path(reactant)
    outfile = str(Path(work_dir) / f"{temp_path.stem}_shifted{temp_path.suffix}")
    interface.pymol.general_cmd(session,[
        ('alter_state', '1', '(all)', f"x+={shift[0]}"), 
        ('alter_state', '1', '(all)', f"y+={shift[1]}"), 
        ('alter_state', '1', '(all)', f"z+={shift[2]}"), 
        ('save', outfile)
    ])

    return outfile

def _place_alphafill(molfile: str,
                     reactant: str,
                     code: str,
                     work_dir: str,
                     similarity_cutoff: float = 0.75,
                     clash_radius: float = 2.0,
                     **kwargs) -> str:
    """Implementation function that places a specified ligand into a molfile structure using the alphafill 
    algorithm (https://doi.org/10.1038/s41592-022-01685-y). SHOULD NOT be called directly by the user. 
    Instead, call `enzy_htp.preparation.place_ligand.place_ligand()`. Assumes that the input ligand has PDB style naming 
    and is in the .mol2 format. When scanning alphafill transplants, the algorithm looks for exact atom label
    matches. If no exact matches exist, then it looks for candidates where the reactant's atom names are a subset 
    of the candidates. If no satisfactory matches are found, it looks for candidates where the candidate's atom
    names are a subset of the reactants. For each subcategory, candidates template molecules are ranked by the 
    number of clashes with the protein structure, with fewer being better. The satisfactory template molecule is
    saved to a .mol2 file in the work_dir and returned. If no template is found, an empty str() is returned.


    Args:
        molfile:
        reactant:
        code:
        work_dir:
        similarity_cutoff:
        clash_radius:

    Returns:
        THe path to the placed reactant in a .mol2 file.

    """
    #TODO(CJ): get the kwargs unpacked here
    _LOGGER.info(f"Beginnning placement of ligand {reactant} into {molfile}...")
    _LOGGER.info(f"Calling out to AlphaFill to fill structure...")
    filled_structure: str = interface.pymol.fetch(code, work_dir)
    filled_structure = interface.alphafill.fill_structure(filled_structure)
    _LOGGER.info("Filled structure using AlphaFill!")

    session = interface.pymol.new_session()
    interface.pymol.general_cmd(session, [("load", filled_structure), ("remove", "solvent"), ("load", molfile),
                                          ("align", Path(filled_structure).stem, Path(molfile).stem)])
    f_df: pd.DataFrame = interface.pymol.collect(session,
                                                 'memory',
                                                 "chain resn resi segi name elem x y z".split(),
                                                 sele=f"{Path(filled_structure).stem} and (not elem H)")
    p_df: pd.DataFrame = interface.pymol.collect(session,
                                                 'memory',
                                                 "chain resn resi segi name elem x y z".split(),
                                                 sele=f"{Path(molfile).stem} and (not elem H)")
    interface.pymol.general_cmd(session, [("delete", Path(molfile).stem)])

    protein_tks = list()
    candidates = list()
    for (chain, resi, resn, segi) in set(zip(f_df.chain, f_df.resi, f_df.resn, f_df.segi)):
        if resn in eh.chemical.THREE_LETTER_AA_MAPPER:
            protein_tks.append(f"( chain {chain} and resi {resi} and resn {resn} and segi {segi})")
            continue

        if resn.upper() in "MG ZN CL NA".split():  #TODO(CJ): make this a function in chemical
            continue
        temp = f_df[(f_df.chain == chain) & (f_df.resi == resi) & (f_df.resn == resn) & (f_df.segi == segi)].reset_index(drop=True)
        candidates.append(deepcopy({'df': temp, 'names': set(temp.name.to_list()), 'key': (chain, resi, resn, segi)}))

    protein_sele: str = "not elem H and " + " or ".join(protein_tks)
    interface.pymol.general_cmd(session, [("delete", "all"), ("load", filled_structure), ("remove", "solvent")])
    #p_df:pd.DataFrame=interface.pymol.collect(session, 'memory', "chain resn resi segi name elem x y z".split(), sele=protein_sele)

    for cc in candidates:
        cc['clashes'] = _count_clashes(p_df, cc['df'], clash_radius)

    interface.pymol.general_cmd(session, [("delete", "all"), ("load", reactant), ("remove", "solvent")])
    r_df: pd.DataFrame = interface.pymol.collect(session, 'memory', "chain resn resi segi name elem x y z".split(), sele="not elem H")

    rct = {'df': r_df, 'names': set(r_df.name.to_list())}

    exact_matches = list()
    rct_is_subset = list()
    can_is_subset = list()
    
    #TODO(CJ): fix this part
    for cc in candidates:
        if rct['names'] == cc['names']:
            exact_matches.append(cc)
            continue

        if rct['names'].issubset(cc['names']):
            rct_is_subset.append(cc)
            continue

        if cc['names'].issubset(rct['names']):
            can_is_subset.append(cc)
            continue

    template = str()
    if exact_matches:
        exact_matches.sort(key=lambda dd: dd['clashes'])
        row = exact_matches[0]['df'].iloc[0]

        _LOGGER.info(f"Found exact match at {row.chain}.{row.resi}.{row.resn}! Transplanting")
        temp_path = Path(reactant)
        outfile = temp_path.parent / f"{temp_path.stem}_placed.mol2"
        template: str = f"{work_dir}/template.mol2"
        interface.pymol.general_cmd(
            session, [("delete", "all"), ("load", filled_structure), ("load", molfile),
                      ("align", Path(filled_structure).stem, Path(molfile).stem),
                      ("save", template,
                       f"{Path(filled_structure).stem} and chain {row.chain} and resn {row.resn} and resi {row.resi} and segi {row.segi}")])
        return align_ligand(template, reactant, outfile=outfile)

    assert False
    if rct_is_subset:
        pass


def _count_clashes(df1: pd.DataFrame, df2: pd.DataFrame, cutoff: float) -> int:
    """Counts the number of clashes between the points in df1 and df2

    Args:
        df1:
        df2:
        cutoff:

    Returns:
        pass
    """
    #TODO(CJ): do checks on the input pandas dataframes
    points1 = np.transpose(np.array([df1.x.to_numpy(), df1.y.to_numpy(), df1.z.to_numpy()]))
    points2 = np.transpose(np.array([df2.x.to_numpy(), df2.y.to_numpy(), df2.z.to_numpy()]))

    count: int = 0

    for pp1 in points1:
        for pp2 in points2:
            count += np.sqrt(np.sum((pp1 - pp2)**2)) <= cutoff

    return count
