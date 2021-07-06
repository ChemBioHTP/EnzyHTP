'''
Misc helper func and class
'''
from AmberMaps import Resi_list
import os
import numpy as np
'''
====
Tree
====
'''
class Child():
    def __init__(self):
        self.parent = None
    def set_parent(self, parent_obj):
        self.parent = parent_obj

        return self

'''
Text
'''
line_feed = '\n'

'''
func
'''
def mkdir(dir):
    if os.path.exists(dir):
        pass
    else:
        os.makedirs(dir)


'''
math
'''
def set_distance(p1,p2,d):
    '''
    return a coord of p3
    -- p3 --
    origin:     p1
    direction: (p1,p2)
    distance:   d
    '''
    p1 = np.array(p1)
    p2 = np.array(p2)
    #direction vector
    v1 = (p2 - p1)/np.linalg.norm(p1 - p2)
    p3 = p1 + v1 * d

    return tuple(p3)

def get_field_strength(p0, c0, p1, p2=None, d1=None):
    '''
    return field strength E of *p0(c0)* at *p1* in direction of *p2-p1* or *d1*
    -- E = kq/r^2 -- (Unit: kcal/mol)
    point charge:   c0 in p0 
    point:          p1
    direction:      p2-p1 or d1
    '''
    # Unit
    k = 332.4   # kcal*Ang/(mol*e^2) = (10^10)*(4.184^-1)*(Na)*(10^-3)*(1.602*10^-19)^2 * 9.0*10^-9 N*m^2/C^2
    q = c0                      # e
    p0 = np.array(p0)           # Ang 
    p1 = np.array(p1)           # Ang
    if d1 == None:
        d1 = np.array(p2) - p1
    else:
        d1 = np.array(d1)
    d1 = d1/np.linalg.norm(d1)  # Ang

    # Get r
    r = p1 - p0
    r_m = np.linalg.norm(r)
    r_u = r/r_m

    # Get E
    E = (k * q / (r_m**2)) * r_u
    # Get E in direction
    Ed = np.dot(E, d1)          # kcal/(mol*e*Ang)

    return Ed


def get_center(p1, p2):
    '''
    return the center of p1 and p2
    '''
    p3 = 0.5 * (np.array(p1) + np.array(p2))

    return tuple(p3)


'''
misc
'''
def decode_atom_mask(stru, mask, ifsolvent=0):
    '''
    decode atom mask and return a list of atom ids.
    Base on the correponding structure obj.
    ===Only support residues in whole now & do not discrimminate chains start id count from 1===
    '''
    atom_ids = []

    resi_ids = mask[1:].strip().split(',')
    for i in range(len(resi_ids)-1, -1, -1):
        if '-' in resi_ids[i]:
            r1 = int(resi_ids[i].split('-')[0])
            r2 = int(resi_ids[i].split('-')[1])
            resi_ids.extend(list(range(r1, r2+1)))
            del resi_ids[i]
    # clean
    resi_ids = [int(i) for i in resi_ids]
    resi_ids.sort()

    all_resi = stru.get_all_residue_unit(ifsolvent=ifsolvent)
    for r_id in resi_ids:
        for resi in all_resi:
            if resi.id == r_id:
                for atom in resi:
                    atom_ids.append(atom.id)

    return atom_ids