import sys
import random

__doc__='''
This module utilize tLEaP to build random mutated structures
------------------------------------------------------------
PDB2Leap()
Input:  PDB file    (The original PDB file // standard Amber format)
        MutaFlag    (The flag of mutation. e.g. E92K)
Output: PDB file    (The structure after mutation // standard Amber format)
------------------------------------------------------------
FlagGen()
Input:  PDB file    (The original PDB file // standard Amber format)
Output: A Random MutaFlag
------------------------------------------------------------
'''

def PDB2Leap(init_PDB_path, MutaFlag):
    
    # Decode the Flag
    Init_resi=MutaFlag[0]
    resi_Index=MutaFlag[1:3]
    Muta_resi=MutaFlag[3]

    #Operate the PDB
    with open(init_PDB_path,'r') as f:
        for line in f:
            if line[25] == resi_Index:
                pass


def FlagGen(init_PDB_path):
    pass

# This part is for test only
PDB1_path=r'C:\Users\shaoqz\OneDrive\Zhongyue\wkFlow\tleap_test\random\2kz2init_amb.pdb'
PDB2Leap(PDB1_path, 'E92K')
# This part is for test only
