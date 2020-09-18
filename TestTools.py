import os

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
Input:  PDB file (after mutation // standard Amber format)
Output: A message of bad contact 
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

def Check_PDB(PDB_path):
    