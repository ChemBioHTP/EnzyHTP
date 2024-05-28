"""Module for fix missing structures of part of the sequence (normally a loop)
or some side-chain heavy atoms.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-09-22
"""
import requests
from typing import List, Union, Dict, Callable

from enzy_htp.chemical import SeqRes

from enzy_htp.structure import Structure

from enzy_htp import interface, config
from enzy_htp.core import _LOGGER
from enzy_htp.core import file_system as fs

def add_missing_residues( stru: Structure, 
                                missing_residues:List[SeqRes], 
                                method:str='modeller', 
                                work_dir:str=None,
                                **kwargs
                                ) -> None:


    """Add in real amino acids with 3D coordinates into a supplied Structure() missing known
    amino acids. This is the main client function that should be called by users. Note that all
    missing Residue()'s are added to the supplied Structure() in place.

    Args:
        stru: The Structure() to which the missing Residue()'s will be added.
        missing_residues: A List[SeqRes] of missing residues, often created by identify_missing_residues().
        method: Which package should be used? 
        work_dir: The directory where the work will take place.

    Returns:
        Nothing.        
    """
    if not work_dir:
        work_dir = config['system.SCRATCH_DIR']

    fs.safe_mkdir( work_dir )

    if not missing_residues:
        return 

    func =  RESIDUE_ADDER_MAPPER.get( method )

    if func is None:
        err_msg:str = f"The supplied method '{method}' is not a valid way to add missing residues. Allowed methods include: {', '.join(RESIDUE_ADDER_MAPPER.keys())}."
        _LOGGER.error( err_msg )
        raise ValueError( err_msg )

    func( stru, missing_residues, work_dir=work_dir, **kwargs )

def identify_missing_residues( code:str ) -> List[SeqRes]:
    """Helper function that gets the missing residues for a given PDB code. Requires
    access to the internet to work. Validation is not directly performed on the code, but will 
    presumably not work with an inavlid one.

    Args:
        code: The 4-letter PDB code to look up.

    Returns:
        A List[SeqRes] of the missing residues in the given PDB code.
    """
    
    url:str=f"https://files.rcsb.org/download/{code}.pdb"
    marker:str='REMARK 465'
    r = requests.get(url)
    if not (r.text.count('ATOM') + r.text.count("HETATM")):
        err_msg:str=f"Unable to get PDB Structure for {code}. Check that this is the right code!"
        _LOGGER.error( err_msg )
        raise ValueError( err_msg )
    
    lines:List[str] = map(lambda ll: ll.decode('utf-8'), r.content.splitlines())
    lines:List[str] = list(filter(lambda ll: ll.startswith(marker), lines))
    lines.reverse()

    if not lines:
        return list()
    
    ll = lines.pop() 
    while ll.find('M RES C SSSEQI') == -1:
        ll = lines.pop()


    result:List[SeqRes] = list()

    for ll in lines:

        three_letter:str=ll[15:18].strip()
        chain_name:str=ll[19].strip()
        idx:int=int(ll[21:26])
        result.append(SeqRes(
            model=1,
            chain=chain_name,
            idx=idx,
            name=three_letter,
            missing=True,
            seq_idx=None))

    result.sort()

    return result

RESIDUE_ADDER_MAPPER:Dict[str, Callable] = {
    'modeller':interface.modeller.add_missing_residues,
    'rosetta':interface.rosetta.add_missing_residues
}
"""TODO(CJ)"""
