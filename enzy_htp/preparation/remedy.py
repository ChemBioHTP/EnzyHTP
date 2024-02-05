"""Module for fix missing structures of part of the sequence (normally a loop)
or some side-chain heavy atoms.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2022-09-22
"""
import requests
from typing import List, Union

from enzy_htp.chemical import SeqRes

from enzy_htp.structure import Structure

from enzy_htp import interface

def add_missing_residues( stru: Structure, 
                                missing_residues:List[SeqRes], 
                                method:str='modeller', 
                                work_dir:str=None,
                                inplace:bool=True
                                ) -> Union[None,Structure]:



    interface.modeller.add_missing_residues( stru, missing_residues, work_dir=work_dir, inplace=inplace)

def identify_missing_residues( code:str ) -> List[SeqRes]:

    #TODO(CJ): add in some checking about the supplkied PDB code. 

    url:str=f"https://files.rcsb.org/download/{code}.pdb"
    marker:str='REMARK 465'
    r = requests.get(url)
    if not r.ok:
        pass
    lines:List[str] = map(lambda ll: ll.decode('utf-8'), r.content.splitlines())
    lines:List[str] = list(filter(lambda ll: ll.startswith(marker), lines))
    lines.reverse()
    
    ll = lines.pop() 
    while ll.find('M RES C SSSEQI') == -1:
        ll = lines.pop()


    result:List[SeqRes] = list()

    for ll in lines:

        three_letter:str=ll[15:18].strip()
        chain_name:str=ll[19].strip()
        idx:int=int(ll[21:26])
        result.append(SeqRes(model=1,chain=chain_name,idx=idx,name=three_letter,missing=True, seq_idx=None))

    result.sort()

    return result
