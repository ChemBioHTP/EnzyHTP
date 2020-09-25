import os,re
import numpy as np
from AmberMaps import *

__doc__='''
This module provides tools to test the whole workflow
------------------------------------------------------------
PDBMin()
Input:  PDB file (standard Amber format)
Output: PDB file (after minimization // standard Amber format)
------------------------------------------------------------
PDB2FF()
Input:  PDB file (standard Amber format)
Output: prmtop & inpcrd 
------------------------------------------------------------
Check_PDB()
Input:  PDB file (after 2FF or minimization // standard Amber format)
Output: A message card of bad contact 
------------------------------------------------------------
'''

def PDB2FF(PDB_path):
    PDB_name=PDB_path[:-4]
    out3_PDB_path=PDB_name+'_water.pdb'
    #make tleap input
    os.system('mkdir tleap_cache')
    tleap_input=open('tleap_ff.in','w')
    tleap_input.write('source leaprc.protein.ff14SB\n')
    tleap_input.write('source leaprc.water.tip3p\n')
    tleap_input.write('a = loadpdb '+PDB_path+'\n')
    tleap_input.write('solvatebox a TIP3PBOX 8\n')
    tleap_input.write('addions a Na+ 0\n') #note here charge should be determined when dealing with more enzyme (mark:enzyme)
    tleap_input.write('saveamberparm a '+PDB_name+'.prmtop '+PDB_name+'.inpcrd\n')
    tleap_input.write('savepdb a '+out3_PDB_path+'\n')
    tleap_input.write('quit\n')
    tleap_input.close()

    #run
    os.system('tleap -s -f tleap_ff.in > tleap_ff_'+PDB_name+'.out')
    os.system('mv *leap_ff* leap.log tleap_cache')

    return (PDB_name+'.prmtop',PDB_name+'.inpcrd',out3_PDB_path)



def PDBMin(PDB_path):
    
    min_cycle='2000'

    PDB_name=PDB_path[:-4]
    out4_PDB_path=PDB_name+'_min.pdb'

    #run pdb2ff
    ff_files=PDB2FF(PDB_path)

    #make sander input
    os.system('mkdir min_cache')
    min_input=open('min.in','w')
    min_input.write('Minimize\n')
    min_input.write(' &cntrl\n')
    min_input.write('  imin=1,\n')
    min_input.write('  ntx=1,\n')
    min_input.write('  irest=0,\n')
    min_input.write('  maxcyc='+min_cycle+',\n')
    min_input.write('  ncyc=10000,\n')
    min_input.write('  ntpr=1000,\n')
    min_input.write('  ntwx=0,\n')
    min_input.write('  cut=8.0,\n')
    min_input.write(' /\n')
    min_input.close()
   

    #run
    os.system('$AMBERHOME/bin/sander -O -i min.in -o min.out -p '+ff_files[0]+' -c '+ff_files[1]+' -r min.rst')

    #rst2pdb
    os.system('ambpdb -p '+ff_files[0]+' -c min.rst > '+out4_PDB_path)
    os.system('mv min.rst min_cache/min_'+PDB_name+'.rst')
    os.system('mv min.out min.in '+ff_files[0]+' '+ff_files[1]+' min_cache')

    return out4_PDB_path

def Check_PDB(PDB_path,MutaFlag):

    Coord=[]
    Target_list=[]
    error_card='Fine\n'

    #decode the Flag
    Init_resi=MutaFlag[0]
    resi_Index=re.search('[0-9]+',MutaFlag).group()
    Muta_resi=MutaFlag[-1]

    # Read coordinate and the target section // super format dependent // Only for minization or 2FF
    with open(PDB_path) as f:
        line_index=1
        Coord_index=0
        for line in f:
            #Skip some line
            skip_flag=0
            #skip the CRYST1 line
            if line_index==1:
                line_index=line_index+1
                continue
            #skip the TER
            if line[:3] == 'TER':
                line_index=line_index+1
                continue
            #skip the END
            if line[:3] == 'END':
                line_index=line_index+1
                continue

            #skip the water and ion(Append in the future)
            skip_list=['Na+','CA','WAT']
            for i in skip_list:
                if line.split()[3] == i:
                    line_index=line_index+1
                    skip_flag=1
                    break
            if skip_flag:
                line_index=line_index+1
                continue

            #Get coordinate form PDB file. The list index is correspond to the (atom index - 1)
            Coord = Coord + [[float(line.split()[5]),float(line.split()[6]),float(line.split()[7])],]

            #Get the target atom list by resi_index
            if line.split()[4] == resi_Index:
                Atom_type=line.split()[2]
                Target_list=Target_list+[(Atom_type,Coord_index),]
            
            Coord_index = Coord_index+1 #index of next element
            line_index=line_index+1 

    # Check each atom in the list by distance
    for atom in Target_list:

        P_A = np.array(Coord[atom[1]])

        #Check each atom expect what is connected
        for i in range(len(Coord)):
            P_B = np.array(Coord[i])
            D_AB = np.linalg.norm(P_B-P_A)
            error_sec=''
            skip_flag=0
            
            #atoms connected (Update to topology version in the future)
            for j in Target_list:                 
                if i == j[1]:
                    # Do not check bond in this version
                    skip_flag=1
            if skip_flag:
                continue
            
            #atoms not connected
            if D_AB <= 0.8:
                error_sec='Bad contact between: atom '+str(i+1)+' and '+str(atom[1]+1)+'\n'

                error_card=error_card+error_sec

    # Check for ring containing Residue
    # How?

    return error_card
    

    

#TestOnly
#print(Check_PDB('2kz2_E92W.rst.pdb','E92W'))
# for i in glob.glob('*min.pdb'):
#     Flag = i.split('.')[0].split('_')[2]
#     print('-------'+Flag+'----------')
#     print(Check_PDB(i,Flag))