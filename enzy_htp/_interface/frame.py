"""Class that encodes the enzy_htp.molecular_mechanics.Frame() class. This datastructure represents
a specific snapshot in time of an enzyme system. Submodule also includes functions to read the Frame()'s
in from .pdb files containing multiple models.


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 06-29-2022
"""
from copy import deepcopy
from typing import List, Dict
from collections import defaultdict

import numpy as np
import pandas as pd

import enzy_htp.chemical as chem
from enzy_htp.core import file_system as fs

class FrameAtom:
    def __init__(self, idx:str, aname:str, resnum:int, resname:str, point, connections:List ):
        self.idx = idx
        self.aname = aname
        self.resnum = resnum
        self.resname = resname 
        self.point = point
        self.connections = connections
        self.new_idx = -1
        self.new_connections = list()
        self.charge = -1

    def update_idx(self, mapper) -> None:
        self.new_idx = mapper[self.idx]
    
    def update_connections(self, mapper) -> None:
        for oc in self.connections:
            self.new_connections.append(mapper[oc])

    def atom_line(self) -> str:
        #TODO(CJ): better way to sanitize the name
        aname = str()
        aname = self.aname[0]
        return f"{aname}  {self.point[0]:.6f} {self.point[1]:.6f} {self.point[2]:.6f}"

def rescale(arr, magnitude):
    arr /= np.linalg.norm(arr)
    return arr*magnitude

class Frame:

    def __init__(self, atoms:List[FrameAtom]):
        self.atoms = atoms

    def update_charges(self, charges):
        assert len(charges) == len(self.atoms)
        for idx, charge in enumerate( charges ):
            self.atoms[idx].charge = charge

    def get_atom_lines(self) -> List[str]:
        result = list()
        for aa in self.atoms:
            result.append( aa.atom_line() )
        return result

    def get_connectivity_lines(self) -> List[str]:
        pass

    def filter_atoms(self, mask) -> Dict:
        assert mask.startswith(':')
        keep_res = list(map(int,mask[1:].split(',')))
        keep, toss = [], dict() 
        for fa in self.atoms:
            if fa.resnum in keep_res:
                keep.append(fa)
            else:
                toss[fa.idx] = fa
        self.atoms = keep
        return toss

    def atoms_needing_caps(self, toss):
        result = dict()        
        for aidx,aa in enumerate(self.atoms):
            for cnt in aa.connections:
                if cnt in toss:
                    result[aidx] = toss[cnt]
        return result

    def apply_H_caps(self, cap_mapper : Dict[int,FrameAtom]):
        cap_atoms = list()
        parent_idx: int
        child : FrameAtom
        for parent_idx,child in cap_mapper.items():
            parent:FrameAtom = self.atoms[parent_idx]
            bond_length = chem.get_h_bond_length(parent.aname)
            if bond_length < 0 :
                raise TypeError()
            start = parent.point
            vec = child.point - parent.point
            vec = rescale(vec, bond_length)
            cap_atoms.append( FrameAtom(child.idx, 'H', child.resnum, 'CAP', start + vec, [parent.idx]))
        self.atoms.extend(cap_atoms)

    def reset_indices(self):

        self.atoms.sort(key=lambda fa: fa.idx)
        idxs = list(map(lambda fa: fa.idx, self.atoms))
        mapper = dict(zip(idxs,range(1,len(idxs)+1)))
        
        atom: FrameAtom
        for atom in self.atoms:
            atom.update_idx( mapper)     
            atom.update_connections( mapper)     

    def apply_mask(self, mask:str, cap:str) -> None:
        toss = self.filter_atoms(mask)
        cap_mapper = self.atoms_needing_caps(toss)
        #TODO(CJ): put in some mask checks and parsing
        if cap is None:
            pass
        elif cap == 'H':
        	self.apply_H_caps(cap_mapper)

        self.reset_indices()

def parse_raw_connections(lines:List[str]) -> Dict:
    """TODO"""
    result = dict()
    
    for ll in lines:
        if not ll.startswith('CONECT'):
            continue
        values = list(map(int,ll.split()[1:]))
        result[values[0]] = values[1:]

    return result


def parse_raw_pdb(lines:List[str]) -> pd.DataFrame:
    """TODO"""
    result = defaultdict(list)
    for ll in lines:
        if not ll.startswith('ATOM') or ll.startswith('HETATM'):
            continue
        result['idx'].append(int(ll[6:11]))
        result['aname'].append(ll[12:16].strip())
        result['resnum'].append(int(ll[22:26]))
        result['resname'].append(ll[17:20])
        result['point'].append(np.array([float(ll[30:38]), float(ll[38:46]), float(ll[46:54])]))
    return pd.DataFrame(result)


def read_charge_list(prmtop:str) -> List[float]:
    lines:List[str] = fs.lines_from_file(prmtop)
    idx = 0
    while not lines[idx].startswith('%FLAG CHARGE'):
        idx+=1
    start = idx + 2
    end = start
    while not lines[end].startswith('%'):
        end += 1
    raw_points = ' '.join(lines[start:end])
    result = list(map(float,raw_points.split()))
    result = [c/18.2223 for c in result]
    return result

def frames_from_pdb(pdb:str) -> List[Frame]:
    """TODO"""
    def desired_line(line:str) -> bool:
        for start in "ATOM HETAM TER END CONECT MODEL".split():
            if line.startswith(start):
                return True
        return False

    def separate(lines:List[str]) -> List[List[str]]:
        result = []
        lines.reverse()
        temp = []
        line = str()
        while len(lines):
            line = lines.pop()      
            if line.startswith('MODEL'):
                continue
            temp.append(line)
            if line.startswith('END'):
                result.append(temp)
                temp = list()

        return (result[:-1], result[-1])

    lines:List[str] = list(filter(desired_line,fs.lines_from_file(pdb)))
    (raw_frames,raw_connects) = separate(lines)
    c_mapper = parse_raw_connections( raw_connects )
    result = list()
    for rf in raw_frames:
        atoms = []
        for i,row in parse_raw_pdb(rf).iterrows():
            temp = deepcopy(row.to_dict())
            temp['connections'] = c_mapper[row.idx]
            atoms.append(FrameAtom(**temp))
        result.append(Frame(atoms))
    return result



#def calculate_electric_fields(frames:List[Frame], 
