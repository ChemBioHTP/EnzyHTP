"""

Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-11-06
"""
import numpy as np

#TODO(CJ): update all this stuff

def calculate_efield( df, p1, direction ):
    result = 0

    k = 332.4

    for i, row in df.iterrows():
        p0 = np.array([row.x, row.y, row.z])
        r = p1 - p0
        r_m:float = np.linalg.norm(r)
        r_u = r/r_m
        result += np.dot((((k*row.charge)/(r_m**2))*r_u), direction)

    return result




def get_field_strength(prmtop:str, frame:str, key1:str, key2:str, location:str='center') -> float: 
   """
   """
   #TODO(CJ): need to check the prmtop file type
   #TODO(CJ): need to check the frame style
   pass
