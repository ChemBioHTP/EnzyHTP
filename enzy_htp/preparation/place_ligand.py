"""Provides functionality for placing a ligand into a protein. All functionalty should be 

Author: Chris Jurich <chris.jurich@vanderbit.edu>
Date: 2023-09-14
"""

from typing import Tuple, Dict


import numpy as np
import pandas as pd
import enzy_htp as eh
from pathlib import Path
from copy import deepcopy
from collections import defaultdict

from enzy_htp import interface, config
from enzy_htp.core import file_system as fs

from enzy_htp.core import _LOGGER

from enzy_htp.structure import PDBParser, Structure

from .align_ligand import align_ligand

def place_ligand(molfile:str, ligand:str, code:str=None, method:str=None, new_res_key:Tuple=None, work_dir:str=None, outfile:str=None, constraints=None, **kwargs ) -> Tuple[str, str]:
    """Algorithm for placing a ligand into  

    Args:
        molfile:
        ligand:
        code:
        method:
        new_res_key:
        work_dir:
        outfile:
        constraints:


    Returns:
        A Tuple[str,str] with the layout (enzyme-ligand complex in .pdb format, exact ligand conformation in .mol2 format).
    """

    if not work_dir:
        work_dir = config['system.SCRATCH_DIR']
    
    fs.safe_mkdir( work_dir )

    if method not in PLACEMENT_MAPPER:
        _LOGGER.error(f"The supplied method '{method}' is not supported. Supported placement methods include: {', '.join(PLACEMENT_MAPPER.keys())}. Exiting...")
        exit( 1 )

    placed_ligand:str=PLACEMENT_MAPPER[method](molfile, ligand, code=code, work_dir=work_dir, constraints=constraints, **kwargs)

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

    args.extend([
        ('save', placed_ligand),
        ('load', molfile),
        ('save', outfile)
    ])
    
    interface.pymol.general_cmd(session, args)

    return (outfile, placed_ligand)




def _place_alphafill(molfile:str, reactant:str, code:str, work_dir:str, similarity_cutoff:float=0.75, clash_radius:float=2.0, **kwargs) -> str:
    """Implementation function that places a specified ligand into a molfile structure using the alphafill 
    algorithm (https://doi.org/10.1038/s41592-022-01685-y). Assumes that the input ligand has PDB style naming 
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
    _LOGGER.info(f"Beginnning placement of ligand {reactant} into {molfile}...")
    _LOGGER.info(f"Calling out to AlphaFill to fill structure...")
    filled_structure:str=interface.pymol.fetch(code, work_dir)
    filled_structure = interface.alphafill.fill_structure(filled_structure)
    _LOGGER.info("Filled structure using AlphaFill!")

    session = interface.pymol.new_session()
    interface.pymol.general_cmd(session, [("load", filled_structure),("remove", "solvent"), ("load", molfile), ("align", Path(filled_structure).stem, Path(molfile).stem)])
    f_df:pd.DataFrame=interface.pymol.collect(session, 'memory', "chain resn resi segi name elem x y z".split(), sele=f"{Path(filled_structure).stem} and (not elem H)")
    p_df:pd.DataFrame=interface.pymol.collect(session, 'memory', "chain resn resi segi name elem x y z".split(), sele=f"{Path(molfile).stem} and (not elem H)")
    interface.pymol.general_cmd(session, [("delete", Path(molfile).stem)])
    
    protein_tks = list()
    candidates = list()
    for (chain, resi, resn, segi) in set(zip(f_df.chain, f_df.resi, f_df.resn, f_df.segi)):
        if resn in eh.chemical.THREE_LETTER_AA_MAPPER:
            protein_tks.append(f"( chain {chain} and resi {resi} and resn {resn} and segi {segi})")
            continue

        if resn.upper() in "MG ZN CL NA".split(): #TODO(CJ): make this a function in chemical
            continue
        temp = f_df[(f_df.chain==chain)&(f_df.resi==resi)&(f_df.resn==resn)&(f_df.segi==segi)].reset_index(drop=True)
        candidates.append(
            deepcopy({'df':temp, 'names':set(temp.name.to_list()), 'key':(chain, resi, resn, segi)})
        )
    
    protein_sele:str = "not elem H and " + " or ".join(protein_tks)
    interface.pymol.general_cmd(session, [("delete", "all"),("load", filled_structure),("remove", "solvent")])
    #p_df:pd.DataFrame=interface.pymol.collect(session, 'memory', "chain resn resi segi name elem x y z".split(), sele=protein_sele)

    for cc in candidates:
        cc['clashes'] = _count_clashes(p_df, cc['df'], clash_radius)

    interface.pymol.general_cmd(session, [("delete", "all"),("load", reactant),("remove", "solvent")])
    r_df:pd.DataFrame=interface.pymol.collect(session, 'memory', "chain resn resi segi name elem x y z".split(), sele="not elem H")

    rct = {'df':r_df, 'names':set(r_df.name.to_list())}

    exact_matches = list()
    rct_is_subset = list()
    can_is_subset = list()

    for cc in candidates:
        if rct['names'] == cc['names']:
            exact_matches.append( cc )            
            continue

        if rct['names'].issubset(cc['names']):
            rct_is_subset.append( cc )
            continue

        if cc['names'].issubset(rct['names']):
            can_is_subset.append( cc )
            continue

    template=str()
    if exact_matches:
        exact_matches.sort(key=lambda dd: dd['clashes'])
        row=exact_matches[0]['df'].iloc[0]
        
        _LOGGER.info(f"Found exact match at {row.chain}.{row.resi}.{row.resn}! Transplanting" )
        temp_path = Path(reactant)
        outfile = temp_path.parent / f"{temp_path.stem}_placed.mol2"
        template:str=f"{work_dir}/template.mol2"
        interface.pymol.general_cmd(session, [("delete", "all"),
                ("load", filled_structure),
                ("load", molfile),
                ("align", Path(filled_structure).stem, Path(molfile).stem),
                ("save", template, f"{Path(filled_structure).stem} and chain {row.chain} and resn {row.resn} and resi {row.resi} and segi {row.segi}")])
        return align_ligand(template, reactant, outfile=outfile)

    assert False
    if rct_is_subset:
        pass



def _count_clashes(df1:pd.DataFrame, df2:pd.DataFrame, cutoff:float) -> int:
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

    count:int = 0

    for pp1 in points1:
        for pp2 in points2:
            count += np.sqrt(np.sum((pp1-pp2)**2)) <= cutoff
    
    return count


PLACEMENT_MAPPER:Dict = {
    'alphafill':_place_alphafill
}
"""TODO(CJ)"""

