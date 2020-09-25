import os
import random
import re
from AmberMaps import *
from TestTools import *

__doc__='''
This module utilize tLEaP to build random mutated structures
------------------------------------------------------------
PDB2PDBwLeap()
Input:  PDB file    (The original PDB file // standard Amber format)
        MutaFlag    (The flag of mutation in Amber PDB index e.g. E92K -> E37K)
Output: PDB file    (The structure after mutation // standard Amber format)
------------------------------------------------------------
FlagGen()
Input:  PDB file    (The original PDB file // standard Amber format)
Output: A Random MutaFlag
------------------------------------------------------------
(List)OldAtoms
A list of atoms kept from the initial residue.
------------------------------------------------------------
'''

def PDB2PDBwLeap(init_PDB_path, MutaFlag):
    
    # Decode the Flag
    Init_resi=MutaFlag[0]
    resi_Index=re.search('[0-9]+',MutaFlag).group()
    Muta_resi=MutaFlag[-1]

    # Operate the PDB
    out_PDB_path1=init_PDB_path[:-4]+'_'+MutaFlag+'_p.pdb'
    out_PDB_path2=init_PDB_path[:-4]+'_'+MutaFlag+'.pdb'

    with open(init_PDB_path,'r') as f:
        with open(out_PDB_path1,'w') as of:
            line_index=1
            for line in f:
                try:
                    if line.split()[4] == resi_Index:
                        #This match method may only work with monomer. Ploymer may need break after first found? (mark:enzyme)
                        #keep if match the OldAtom list
                        for i in OldAtoms:
                            if i == line.split()[2]:
                                #change the traget residue name !!WARNING!! strong format limitation here!
                                
                                new_line=line[:17]+Resi_map[Muta_resi]+line[20:]

                                of.write(new_line)
                                break
                        #Do not keep if not match & in the target residue.
                    else:
                        of.write(line)
                
                except IndexError:
                    of.write(line)
                    print('Warning: Not a data line -> '+line+'---'+str(line_index))

                line_index=line_index+1
    
    # Run tLeap
    #make input
    os.system('mkdir tleap_cache')
    leap_input=open('tleap.in','w')
    leap_input.write('source leaprc.protein.ff14SB\n')
    leap_input.write('a = loadpdb '+out_PDB_path1+'\n')
    leap_input.write('savepdb a '+out_PDB_path2+'\n')
    leap_input.write('quit\n')
    leap_input.close()
    #run
    os.system('tleap -s -f tleap.in > tleap.out')
    os.system('mv *leap.* tleap_cache')

    return out_PDB_path2



def FlagGen(init_PDB_path):
    out_flag=''
    tot_resi=0
    resi_1=''
    resi_2=''
    Muta_idx=''

    #This method to count the totle residue in a chain is order depending. The target chain has to be the first.(mark:enzyme)
    with open(init_PDB_path,'r') as f:
        lines=f.readlines()
        for i in range(len(lines)):
            if lines[i].strip() == 'TER':
                tot_resi=int(lines[i-1].split()[4])
                break
        #Generate the mutation index
        Muta_idx=str(random.randint(1,tot_resi))
        #obtain resi_1
        for line in lines:
            if line.split()[4] == Muta_idx:
                resi_1_p=line.split()[3]
                resi_1=Resi_map2[resi_1_p]
                break
        #Generate resi_2
        resi_2=Resi_list[random.randint(0,len(Resi_list)-1)]
        # Check if the same resi
        while resi_2 == resi_1:
            resi_2=Resi_list[random.randint(0,len(Resi_list)-1)]
            


    out_flag=resi_1+Muta_idx+resi_2
    return out_flag

OldAtoms=['N','H','CA','HA','CB','C','O']




# This part is for test only
#PDB1_path='2kz2init_amb.pdb'
#print(FlagGen(PDB1_path))
#PDB2_path=PDB2PDBwLeap(PDB1_path, 'E3K')
#PDBMin(PDB2_path)
# This part is for test only
#TestOnly