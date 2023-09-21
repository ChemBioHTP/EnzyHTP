"""
"""
from pathlib import Path
import pandas as pd

from enzy_htp import interface, config

import enzy_htp.chemical as chem

from enzy_htp.core import file_system as fs


def create_cluster(fname: str, sele_str: str, outfile: str = None, cap_strategy: str = 'H', work_dir: str = None, session=None) -> str:
    """ 
    """

    if work_dir is None:
        work_dir = config['system.WORK_DIR']

    if session is None:
        session = interface.pymol.new_session()

    fs.check_file_exists(fname)

    if outfile is None:
        temp_path = Path(fname)
        outfile: str = f"{work_dir}/{temp_path.stem}_cluster{temp_path.suffix}"

    #TODO(CJ): add file check for fname
    obj_name: str = '__eh_cluster'
    interface.pymol.general_cmd(session, [('load', fname), ('select', sele_str), ('create', obj_name, sele_str)])

    df: pd.DataFrame = interface.pymol.collect(session, 'memory', "chain resi resn name".split())

    args = list()
    for i, row in df.iterrows():
        if row['name'] not in "N C".split():
            continue

        if not row.resn in chem.THREE_LETTER_AA_MAPPER:
            continue

        args.append(('h_add', f"chain {row.chain} and resi {row.resi} and resn {row.resn} and name {row['name']}"))

    args.append(("save", outfile, obj_name))

    interface.pymol.general_cmd(session, args)

    if cap_strategy == 'H':
        return outfile

    if cap_strategy == 'CH3':

        args = []
        orig = set(zip(df.chain, df.resi, df.resn, df['name']))
        updated: pd.DataFrame = interface.pymol.collect(session, 'memory', "chain resi resn name".split())
        new = set(zip(updated.chain, updated.resi, updated.resn, updated['name']))

        for tp in new:
            if tp in orig:
                continue
            if tp[3] == 'H01':
                new_name = 'C21'
            elif tp[3] == 'H02':
                new_name = 'C22'
            else:
                assert False
            args.extend([
                ('alter', f"chain {tp[0]} and resi {tp[1]} and resn {tp[2]} and name {tp[3]}", "elem='C'"),
                ('alter', f"chain {tp[0]} and resi {tp[1]} and resn {tp[2]} and name {tp[3]}", f"name='{new_name}'"),
            ])

        args.extend([('valence', 'guess', 'name C21 or name C22'), ('h_add', 'name C21'), ('h_add', 'name C22'),
                     ("save", outfile, obj_name)])
        interface.pymol.general_cmd(session, args)

    return outfile
