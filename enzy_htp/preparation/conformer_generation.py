"""Contains code for generating ligand conformers. Uses the Ligand() object and adds conformer coordinates directly to the 
object, in-place. Users should access functionality through the generate_conformers() top-level function. It allows users
to use multiple conformer generation softwares including:
    
    + bcl
    + rdkit

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-10-05
"""
from typing import Dict

import numpy as np
import pandas as pd
from rdkit.Chem import AllChem

from enzy_htp.core import file_system as fs
from enzy_htp import interface, config
from enzy_htp.structure import PDBParser, Mol2Parser, Structure, Ligand




def generate_conformers( ligand : Ligand, n_conformers:int=100, engine:str='bcl', work_dir:str=None, **kwargs):
    """Create conformers for the specified ligand and add them to the ligand object. Conformers can be created using multiple
    engines. All conformer addition is done in place to the Ligand and all temp files are deleted. No sanitation of the 
    input Ligand() is performed. A warning is logged if enough conformers are not produced.
    
    Args:
        ligand: The Ligand() object to add conformers to.
        n_conformers: How many conformers to produce. Must be positive and defaults to 100.
        engine: A str() specifying which conformer engine to use. Defaults to 'bcl'.
        work_dir: The directory to store temporary files. Defaults to system.SCRATCH_DIR if not supplied.

    Returns:
        Nothing.
    """
    if n_conformers <= 0:
        _LOGGER.error(f"Supplied n_conformers value ({n_conformers}) is invalid. Must be positive. Exiting...")
        exit( 1 )

    if engine not in METHOD_MAPPER:
        _LOGGER.error(f"The supplied method {engine} is not supported. Allowed engines are {', '.join(METHOD_MAPPER.values())}. Exiting...")
        exit( 1 )

    if work_dir is None:
        work_dir = config['system.SCRATCH_DIR']

    orig_count:int=ligand.n_conformers()
    METHOD_MAPPER[engine](ligand, n_conformers, work_dir)
    updated_count:int=ligand.n_conformers()
    
    n_created:int = updated_count - orig_count

    if n_created < n_conformers:
        _LOGGER.warning(f"Could not make specified number of conformers. Tried to make {n_conformers}, but made {n_created}")
    

def _generate_conformers_bcl( ligand : Ligand, n_conformers:int, work_dir:str, **kwargs) -> None:
    """Implementation method to generate conformers with BCL. SHOULD NOT be called by the user. Instead, call 
    generate_conformers().
    
    Args:
        ligand: The Ligand() object to add conformers to.
        n_conformers: How many conformers to produce.
        work_dir: Directory to store temporary files.

    Returns:
        Nothing.
    """
    to_delete:List[str] = list()

    temp_file:str = f"{work_dir}/__temp.mol2"
    to_delete.append( temp_file )

    fs.safe_rm( temp_file )

    parser = Mol2Parser()
    parser.save_ligand(temp_file, ligand)

    session = interface.pymol.new_session()

    interface.pymol.general_cmd(session,[('load', temp_file)])

    temp_sdf:str=interface.pymol.convert(session, temp_file, new_ext='.sdf')
    conformers:str=interface.bcl.generate_conformers(temp_sdf, n_conformers)
    to_delete.extend([ temp_sdf, conformers])
    
    n_states:int=interface.pymol.general_cmd(session, [('delete', 'all'), ('load', conformers), ('count_states',)])[-1]

    for idx in range( n_states ):
        df:pd.DataFrame=interface.pymol.collect(session, 'memory', "x y z".split(), state=(idx+1))
        locations = list()
        for i, row in df.iterrows():
            locations.append((row.x, row.y, row.z))
        ligand.add_conformer( locations )

    for td in to_delete:
        fs.safe_rm( td )


def _generate_conformers_rdkit( ligand : Ligand, n_conformers:int, work_dir:str, **kwargs ) -> None:
    """Implementation method to generate conformers with rdkit. SHOULD NOT be called by the user. Instead, call 
    generate_conformers().
    
    Args:
        ligand: The Ligand() object to add conformers to.
        n_conformers: How many conformers to produce.
        work_dir: Directory to store temporary files.

    Returns:
        Nothing.
    """
    to_delete:List[str] = list()

    temp_file:str = f"{work_dir}/__temp.mol2"
    to_delete.append( temp_file )

    fs.safe_rm( temp_file )

    parser = Mol2Parser()
    parser.save_ligand(temp_file, ligand)

    mol = interface.rdkit._load_molecule( temp_file )

    confs = AllChem.EmbedMultipleConfs(mol, useRandomCoords=True, numConfs=n_conformers, numThreads=0)

    for conf_idx in list(confs):
        locations = list()
        conformer = mol.GetConformer(conf_idx)
        for aidx, atom in enumerate( mol.GetAtoms() ):
            pos=conformer.GetAtomPosition(aidx)
            locations.append((pos.x, pos.y, pos.z))
        ligand.add_conformer( locations )

    for td in to_delete:
        fs.safe_rm( td )


METHOD_MAPPER:Dict={
    'bcl':_generate_conformers_bcl,
    'rdkit':_generate_conformers_rdkit
}
"""Mapper for the methods available for conformer generation. Implementation detail not to be touched by users."""
