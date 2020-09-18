import os
import re
from random import randint
from AmberMaps import *

__doc__='''
This module defines the potential object to operate. Included some common I/O and methods.
-------------------------------------------------------------------------------------
Class PDB
-------------------------------------------------------------------------------------
self.__init__()
Input:  PDB path (any format)
Output: self.path 
        self.name 
        self.ifformat (Judge by the first line)
-------------------------------------------------------------------------------------
self.get_tot_resi()
Input:  self (standard Amber format)
Output: self.tot_resi
-------------------------------------------------------------------------------------
'''

class PDB(object):

    path=''
    name=''
    ifformat=0

    tot_resi=0

    topology=''
    

    #MutaFlag=''


    #regex pattern of standard Amber format line
    ATOM_pattern=r'ATOM  [ ,0-9][ ,0-9][ ,0-9][ ,0-9][0-9]  [A-Z][ ,0-9,A-Z][ ,0-9,A-Z] [A-Z][A-Z][A-Z]  [ ,0-9][ ,0-9][ ,0-9][0-9]    [ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9][ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9][ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9](?:         [ ,A-Z][ ,A-Z][A-Z])?'
    HETATM_pattern=r'HETATM[ ,0-9][ ,0-9][ ,0-9][ ,0-9][0-9]  [A-Z][ ,0-9,A-Z][ ,0-9,A-Z] [A-Z][A-Z][A-Z]  [ ,0-9][ ,0-9][ ,0-9][0-9]    [ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9][ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9][ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9](?:         [ ,A-Z][ ,A-Z][A-Z])?'

    def __init__(self,PDB_PATH):
        self.path=PDB_PATH
        self.name=self.path.split(os.sep)[-1][:-4]

        #self.ifformat
        #Judge if this is the standard Amber format by just read two line.
        with open(self.path) as f:
            line_index=1
            for line in f:
                if line_index==1:
                    line_index=line_index+1
                    continue
                # check by different length. Need improvement in the future for other target type.                 
                ATOM_if = re.match(self.ATOM_pattern,line) != None
                HATATM_if = re.match(self.HETATM_pattern,line) != None

                self.ifformat = ATOM_if or HATATM_if
                break
        
        if not self.ifformat:
            print('WARNING: This file is not in standard Amber format')

    def get_tot_resi(self):
        #This method of counting the totle residue in a chain is order depending. The target chain has to be the first.(mark:enzyme)
        with open(self.path,'r') as f:
            lines=f.readlines()
            for i in range(len(lines)):
                if lines[i].strip() == 'TER':
                    self.tot_resi=int(lines[i-1].split()[4])
                    break

    


    # def Costum_MutaFlag(self,Flag):
    #     self.MutaFlag=Flag

    # def Random_MutaFlag(self):
    #     resi_1=''
    #     resi_2=''
    #     Muta_idx=''

    #     #This method to count the totle residue in a chain is order depending. The target chain has to be the first.(mark:enzyme)
    #     with open(self.path,'r') as f:
    #         lines=f.readlines()
    #         for i in range(len(lines)):
    #             if lines[i].strip() == 'TER':
    #                 self.tot_resi=int(lines[i-1].split()[4])
    #                 break
    #         #Generate the mutation index
    #         Muta_idx=str(randint(1,self.tot_resi))
    #         #obtain resi_1
    #         for line in lines:
    #             if line.split()[4] == Muta_idx:
    #                 resi_1_p=line.split()[3]
    #                 resi_1=Resi_map2[resi_1_p]
    #                 break
    #         #Generate resi_2
    #         resi_2=Resi_list[randint(0,len(Resi_list)-1)]
    #         # Check if the same resi
    #         while resi_2 == resi_1:
    #             resi_2=Resi_list[randint(0,len(Resi_list)-1)]
                
    #     self.MutaFlag=resi_1+Muta_idx+resi_2



    



#TestOnly
a=PDB(r'C:\Users\shaoqz\OneDrive\Zhongyue\wkFlow\tleap_test\random\2kz2init_amb_A33T_min.pdb')