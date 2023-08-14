"""
"""
import pandas as pd

from enzy_htp import interface, config

import enzy_htp.chemical as chem




def create_cluster( fname:str, outfile:str, sele_str:str, cap_strategy:str='H', work_dir:str=None, session=None) -> str:
    """ 
    """

    if work_dir is None:
        work_dir = config['system.WORK_DIR']

    if session is None:
        session = interface.pymol.new_session()

    #TODO(CJ): add file check for fname
    obj_name:str='__eh_cluster'
    interface.pymol.general_cmd(session, [('load', fname), ('select', sele_str),('create',obj_name,sele_str)])


    df:pd.DataFrame=interface.pymol.collect(session, 'memory', "chain resi resn name".split())


    args = list()
    for i,row in df.iterrows():
        if row['name'] not in "N C".split():
            continue
        
        if not row.resn in chem.THREE_LETTER_AA_MAPPER:
            continue
        
        args.append(
            ('h_add', f"chain {row.chain} and resi {row.resi} and resn {row.resn} and name {row['name']}")
        )

    interface.pymol.general_cmd(session, args )

    outfile = f"{work_dir}/try.pdb"
    interface.pymol.general_cmd(session, [("save", outfile, obj_name)])

    if cap_strategy == 'H':
        return outfile


    if cap_strategy == 'CH3':
        
        args = []
        orig = set(
            zip(df.chain,df.resi,df.resn,df['name'])
        )
        updated:pd.DataFrame=interface.pymol.collect(session, 'memory', "chain resi resn name".split())
        new = set(
            zip(updated.chain,updated.resi,updated.resn,updated['name'])
        )

        for tp in new:
            if tp in orig:
                continue
            if tp[3] == 'H01':
                new_name='C21'
            elif tp[3] == 'H02':
                new_name='C22'
            else:
                assert False
            args.extend([
                ('alter', f"chain {tp[0]} and resi {tp[1]} and resn {tp[2]} and name {tp[3]}", "elem='C'"),
                ('alter', f"chain {tp[0]} and resi {tp[1]} and resn {tp[2]} and name {tp[3]}", f"name='{new_name}'"),
            ])

        outfile = f"{work_dir}/try.mol2"
        args.extend([
            ('valence', 'guess', 'name C21 or name C22'),
            ('h_add','name C21'),
            ('h_add','name C22'),
            ("save", outfile, obj_name)
        ])
        interface.pymol.general_cmd(session, args )


