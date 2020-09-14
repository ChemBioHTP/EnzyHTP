import os
import random
from AmberMaps import *
from TestTools import *#TestOnly

__doc__='''
This module utilize tLEaP to build random mutated structures
------------------------------------------------------------
PDB2Leap()
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

def PDB2Leap(init_PDB_path, MutaFlag):
    
    # Decode the Flag
    Init_resi=MutaFlag[0]
    resi_Index=MutaFlag[1:3]
    Muta_resi=MutaFlag[3]

    # Operate the PDB
    out_PDB_path=init_PDB_path[:-4]+'_'+MutaFlag+'_p.pdb'
    out_PDB_path2=init_PDB_path[:-4]+'_'+MutaFlag+'.pdb'

    with open(init_PDB_path,'r') as f:
        with open(out_PDB_path,'w') as of:
            line_index=1
            for line in f:
                try:
                    if line.split()[4] == resi_Index:
                        #keep if match the OldAtom list
                        for i in OldAtoms:
                            if i == line.split()[2]:
                                #change the traget residue name !!WARNING!! format limitation here!
                                #Another method but still depend on format
                                #Ele=line.split()                            
                                #new_line='{:<8}'.format(Ele[0])+
                                
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
    leap_input.write('a = loadpdb '+out_PDB_path+'\n')
    leap_input.write('savepdb a '+out_PDB_path2+'\n')
    leap_input.write('quit\n')
    leap_input.close()
    #run
    os.system('tleap -s -f tleap.in > tleap.out')
    os.system('mv *leap.* tleap_cache')

    return out_PDB_path2



def FlagGen(init_PDB_path):
    pass

OldAtoms=['N','H','CA','HA','CB','C','O']




# This part is for test only
PDB1_path='2kz2init_amb.pdb'
PDB2_path=PDB2Leap(PDB1_path, 'E37K')
PDBMin(PDB2_path)
# This part is for test only
#TestOnly