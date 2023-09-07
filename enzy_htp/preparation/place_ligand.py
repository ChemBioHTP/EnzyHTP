import numpy as np
import pandas as pd
import enzy_htp as eh
from pathlib import Path
from copy import deepcopy
from collections import defaultdict

from enzy_htp import interface, config

from enzy_htp.structure import PDBParser, Structure


def _map_hydrogens(session):
    
    df:pd.DataFrame = interface.pymol.collect(session, "memory", "name elem x y z".split())
    h_df = df[df.elem=='H'].reset_index(drop=True)
    ha_df = df[df.elem!='H'].reset_index(drop=True)

    ha_names = ha_df['name'].to_list()
    ha_pts = np.transpose(np.array([ha_df.x.to_numpy(), ha_df.y.to_numpy(), ha_df.z.to_numpy()]))
    
    mapper = defaultdict(list)

    for i, row in h_df.iterrows():
        point = np.array([row.x, row.y, row.z])
        
        idx = np.argmin(np.sqrt(np.sum((ha_pts-point)**2,axis=1)))
        mapper[ha_names[idx]].append( row['name'] )
    
    return mapper



def _shift_reactant(outfile:str, reactant:str, df:pd.DataFrame):

    session = interface.pymol.new_session()
    interface.pymol.general_cmd(session, [('load', reactant)])

    orig_h_mapper = _map_hydrogens(session)
    
    args = list()
    for i, row in df.iterrows():
        args.append(
            ('alter_state', -1, f'name {row["name"]}', f'(x,y,z) = ({row["x"]},{row["y"]},{row["z"]})'),
        )
    
    interface.pymol.general_cmd(session, args)
    for heavy_atom, h_names in orig_h_mapper.items():
        args = list()
        for hh in h_names:
            args.append(('remove', f"name {hh}"))
        args.append(('h_add', f"name {heavy_atom}"))
        
        interface.pymol.general_cmd(session, args)
        updated_h_mapper = _map_hydrogens(session)
        args = list()
        for orig, new in zip(sorted(h_names),sorted(updated_h_mapper[heavy_atom])):
            args.append(('alter', f'name {new}', f'name="{orig}"'))
            
        interface.pymol.general_cmd(session, args)


    interface.pymol.general_cmd(session, [('save', outfile)])
    return outfile


def _count_clashes(df1:pd.DataFrame, df2:pd.DataFrame, cutoff=2.0) -> int:
    
    points1 = np.transpose(np.array([df1.x.to_numpy(), df1.y.to_numpy(), df1.z.to_numpy()]))
    points2 = np.transpose(np.array([df2.x.to_numpy(), df2.y.to_numpy(), df2.z.to_numpy()]))

    count = 0

    for pp1 in points1:
        for pp2 in points2:
            count += np.sqrt(np.sum((pp1-pp2)**2)) <= cutoff
    
    return count

def place_reactant(filled_structure:str, reactant:str) -> str:

    session = interface.pymol.new_session()
    interface.pymol.general_cmd(session, [("load", filled_structure),("remove", "solvent")])
    f_df:pd.DataFrame=interface.pymol.collect(session, 'memory', "chain resn resi segi name elem x y z".split(), sele="not elem H")
    
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
    p_df:pd.DataFrame=interface.pymol.collect(session, 'memory', "chain resn resi segi name elem x y z".split(), sele=protein_sele)

    
    for cc in candidates:
        cc['clashes'] = _count_clashes(p_df, cc['df'])

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


    if exact_matches:
        exact_matches.sort(key=lambda dd: dd['clashes'])
        temp_path = Path(reactant)
        outfile = temp_path.parent / f"{temp_path.stem}_placed.mol2"
        return _shift_reactant(outfile, reactant, exact_matches[0]['df'])
            

    assert False
    if rct_is_subset:
        pass

