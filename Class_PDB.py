import os
import re
from random import randint
from AmberMaps import *

__doc__='''
This module defines the potential object to operate. Included some common I/O and methods.
-------------------------------------------------------------------------------------
Class PDB
-------------------------------------------------------------------------------------
__init__(self,path,Flag)
    Input:  PDB path (any format), MutaFlag
    Output: self.path 
            self.name 
            self.ifformat (Judge by the first line)
            self.MutaFlag
          # Connect with the database // add a flow for database?
-------------------------------------------------------------------------------------
Information collecting methods:
-------------------------------------------------------------------------------------
get_tot_resi(self)
    Input:  self (standard Amber format)
    Output: self.tot_resi
get_index_resi(self, index)
get_coord(self, index)
-------------------------------------------------------------------------------------
PDB operating methods: (changes the self.path to indicated new pdb)
-------------------------------------------------------------------------------------
PDB2PDBwLeap(self):
PDBMin(self,cycle):
rm_wat(self): remove water and ion for current pdb. (For potential docking)
-------------------------------------------------------------------------------------
Other Mutation Tools:
-------------------------------------------------------------------------------------
FlagGen(self):
PDB_check(self):
-------------------------------------------------------------------------------------
Input file generating methods:
PDB2FF(self):
PDBMD(self):
-------------------------------------------------------------------------------------
'''

class PDB(object):

    path=''
    name=''
    prmtop_path=''
    inpcrd_path=''
    stage=0 #For debug indicate which stage is the program in
    ifformat=0

    MutaFlag=[]
    tot_resi=0
    Coord=[]

    #current operating residue
    current_index=0
    current_resi_name_3=''
    current_resi_name=''
    current_resi_atom_list=[]
    
    


    #regex pattern of standard Amber format line
    ATOM_pattern=r'ATOM  [ ,0-9][ ,0-9][ ,0-9][ ,0-9][0-9]  [A-Z][ ,0-9,A-Z][ ,0-9,A-Z] [A-Z][A-Z][A-Z]  [ ,0-9][ ,0-9][ ,0-9][0-9]    [ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9][ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9][ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9](?:         [ ,A-Z][ ,A-Z][A-Z])?'
    HETATM_pattern=r'HETATM[ ,0-9][ ,0-9][ ,0-9][ ,0-9][0-9]  [A-Z][ ,0-9,A-Z][ ,0-9,A-Z] [A-Z][A-Z][A-Z]  [ ,0-9][ ,0-9][ ,0-9][0-9]    [ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9][ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9][ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9](?:         [ ,A-Z][ ,A-Z][A-Z])?'


    def __init__(self,PDB_PATH,Flag=''):

        self.path=PDB_PATH
        self.name=self.path.split(os.sep)[-1][:-4]
        if len(Flag) != 0:
            self.MutaFlag = [Flag[0],re.search('[0-9]+',Flag).group(),Flag[-1]]

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
        #This method takes the first chain of a PDB file to count the total number of residues. (mark:enzyme)
        with open(self.path,'r') as f:
            lines=f.readlines()
            for i in range(len(lines)):
                if lines[i].strip() == 'TER':
                    self.tot_resi=int(lines[i-1].split()[4])
                    break

    def get_index_resi(self,index):

        self.current_index = index

        with open(self.path) as f:
            i=0
            for line in f:
                if line.split()[4] == self.current_index:
                    if i==0:
                        self.current_resi_name_3=line.split()[3]
                        self.current_resi_name=Resi_map2[self.current_resi_name_3]
                    Atom_type=line.split()[2]
                    Atom_index=line.split()[1]
                    self.current_resi_atom_list=self.current_resi_atom_list+[(Atom_type,Atom_index),]
                    i=i+1
    
    def get_coord(self):
        # The method to extract coordinate depends critically on the file format. (Only takes PDB after minimization and '2FF')
        self.Coord=[]
        with open(self.path) as f:
            line_index=1
            for line in f:
                #Skip some lines
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

                #Get coordinate form PDB file. The list index is correspond to the (atom index - 1) // may be there will be some inconsistancy? (code:enzyme)
                self.Coord = self.Coord + [[float(line.split()[5]),float(line.split()[6]),float(line.split()[7])],]
                
                line_index=line_index+1 





    def PDB2PDBwLeap(self):
        '''
        Use self.MutaFlag and tLeap to build mutated structure.
        Save changes to self.path and self.name
        self.stage=1
        '''
    
        resi_Index=self.MutaFlag[1]
        Muta_resi=self.MutaFlag[2]
        Flag_name=self.MutaFlag[0]+self.MutaFlag[1]+self.MutaFlag[2]
        OldAtoms=['N','H','CA','HA','CB','C','O']
        #fix for Gly
        if Muta_resi == 'G':
            OldAtoms=['N','H','CA','C','O']

        # Operate the PDB
        out_PDB_path1=self.path[:-4]+'_'+Flag_name+'_tmp.pdb'
        out_PDB_path2=self.path[:-4]+'_'+Flag_name+'.pdb'

        with open(self.path,'r') as f:
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

        #Update the file
        self.path = out_PDB_path2
        self.name = self.path.split(os.sep)[-1][:-4]
        self.stage = 1

        return out_PDB_path2


    def Random_MutaFlag(self):
        '''
        Use self.path & self.tot_resi
        Save changes to self.MutaFlag 
        '''
        resi_1=''
        resi_2=''
        Muta_idx=''

        #This method takes the first chain of a PDB file to count the total number of residues. (mark:enzyme)
        with open(self.path,'r') as f:
            lines=f.readlines()
            for i in range(len(lines)):
                if lines[i].strip() == 'TER':
                    self.tot_resi=int(lines[i-1].split()[4])
                    break
            #Generate the mutation index
            Muta_idx=str(randint(1,self.tot_resi))
            #obtain resi_1
            for line in lines:
                if line.split()[4] == Muta_idx:
                    resi_1_p=line.split()[3]
                    resi_1=Resi_map2[resi_1_p]
                    break
        #Generate resi_2
        resi_2=Resi_list[randint(0,len(Resi_list)-1)]
        # Check if the same resi
        while resi_2 == resi_1:
            resi_2=Resi_list[randint(0,len(Resi_list)-1)]
                
        self.MutaFlag=[resi_1,Muta_idx,resi_2]

    def PDB2FF(self):
        #out3_PDB_path=self.name+'_water.pdb'
        #make tleap input
        os.system('mkdir tleap_cache')
        tleap_input=open('tleap_ff.in','w')
        tleap_input.write('source leaprc.protein.ff14SB\n')
        tleap_input.write('source leaprc.water.tip3p\n')
        tleap_input.write('a = loadpdb '+self.path+'\n')
        tleap_input.write('solvatebox a TIP3PBOX 10\n')
        tleap_input.write('addions a Na+ 0\n')
        tleap_input.write('addions a Cl- 0\n')
        tleap_input.write('saveamberparm a '+self.name+'.prmtop '+self.name+'.inpcrd\n')
        #tleap_input.write('savepdb a '+out3_PDB_path+'\n')
        tleap_input.write('quit\n')
        tleap_input.close()

        #run
        os.system('tleap -s -f tleap_ff.in > tleap_ff_'+self.name+'.out')
        os.system('mv *leap_ff* leap.log tleap_cache')

        self.prmtop_path=self.name+'.prmtop'
        self.inpcrd_path=self.name+'.inpcrd'
        self.stage = 2

        return (self.prmtop_path,self.inpcrd_path)



    def PDBMin(self,cycle='2000'):
        
        out4_PDB_path=self.path[:-4]+'_min.pdb'

        #make sander input
        os.system('mkdir min_cache')
        min_input=open('min.in','w')
        min_input.write('Minimize\n')
        min_input.write(' &cntrl\n')
        min_input.write('  imin=1,\n')
        min_input.write('  ntx=1,\n')
        min_input.write('  irest=0,\n')
        min_input.write('  maxcyc='+cycle+',\n')
        min_input.write('  ncyc=10000,\n')
        min_input.write('  ntpr=1000,\n')
        min_input.write('  ntwx=0,\n')
        min_input.write('  cut=8.0,\n')
        min_input.write(' /\n')
        min_input.close()
    

        #run
        os.system('mpirun -np 8 $AMBERHOME/bin/sander.MPI -O -i min.in -o min.out -p '+self.prmtop_path+' -c '+self.inpcrd_path+' -r min.rst')

        #rst2pdb
        os.system('ambpdb -p '+self.prmtop_path+' -c min.rst > '+out4_PDB_path)
        os.system('mv min.rst min_cache/min_'+self.name+'.rst')
        os.system('mv min.out min.in '+self.prmtop_path+' '+self.inpcrd_path+' min_cache')

        self.path = out4_PDB_path
        self.name = self.name+'_min'
        self.prmtop_path='min_cache/'+self.prmtop_path
        self.inpcrd_path='min_cache/'+self.inpcrd_path
        self.stage = 3

        return out4_PDB_path

    def rm_wat(self):
        '''
        Remove water and ion for the pdb. Remians the same if there's no water or ion.
        Now only skip [Na+,Cl-,WAT] // Append more in the future.
        Save changed files into self.path.
        '''
        out_path = self.path[:-4]+'_rmW.pdb'

        with open(self.path) as f:
            with open(out_path,'w') as of:
                skip_list=['Na+','Cl-','WAT']
                change_flag=0
                for line in f:
                    #Skip some lines
                    skip_flag=0
                    #skip the CRYST1 line
                    if line[:6] == 'CRYST1':
                        change_flag=1
                        continue
                    #keep TER and END
                    if line[:3] == 'TER' or line[:3] == 'END':
                        of.write(line)
                        continue

                    #skip the water and ion(Append in the future)
                    for i in skip_list:
                        if line.split()[3] == i:
                            skip_flag=1
                            break
                    if skip_flag:
                        change_flag=1
                        continue

                    of.write(line)
                if not change_flag:
                    print('No change.')

        self.path=out_path
        self.name=self.name+'_rmW'

        return self.path

    def PDBMD(self,conf_path,tag=''):
        '''
        Use self.prmtop_path and self.inpcrd_path to initilize a MD simulation.
        The MD configuration files are assigned by the conf_path. Be sure to use the same tag when files are taged.
        Return the rst path of the prod step.
        '''
        Flag_name=self.MutaFlag[0]+self.MutaFlag[1]+self.MutaFlag[2]
        out_path=conf_path

        #run sander
        os.system('mpirun -np 8 $AMBERHOME/bin/sander.MPI -O -i '+conf_path+'/min'+tag+'.in -o ' +out_path+'/min_' +Flag_name+'.out -p '+self.prmtop_path+' -c '+self.inpcrd_path+' -r '+out_path+'/min_'+Flag_name+'.rst')
        os.system('mpirun -np 8 $AMBERHOME/bin/sander.MPI -O -i '+conf_path+'/heat'+tag+'.in -o '+out_path+'/heat_'+Flag_name+'.out -p '+self.prmtop_path+' -c '+out_path+'/min_'+Flag_name+'.rst -ref ' +out_path+'/min_'+Flag_name+'.rst -r ' +out_path+'/heat_'+Flag_name+'.rst')
        os.system('mpirun -np 8 $AMBERHOME/bin/sander.MPI -O -i '+conf_path+'/equi'+tag+'.in -o '+out_path+'/equi_'+Flag_name+'.out -p '+self.prmtop_path+' -c '+out_path+'/heat_'+Flag_name+'.rst -ref '+out_path+'/heat_'+Flag_name+'.rst -r '+out_path+'/equi_'+Flag_name+'.rst')
        os.system('mpirun -np 8 $AMBERHOME/bin/sander.MPI -O -i '+conf_path+'/prod'+tag+'.in -o '+out_path+'/prod_'+Flag_name+'.out -p '+self.prmtop_path+' -c '+out_path+'/equi_'+Flag_name+'.rst -ref '+out_path+'/equi_'+Flag_name+'.rst -r '+out_path+'/prod_'+Flag_name+'.rst -x '+out_path+'/prod'+Flag_name+'.nc')

        #Or add a class constant?
        return 'MD/prod_'+Flag_name+'.rst'

    



#TestOnly
# a=PDB(r'2kz2init_amb.pdb','E37K')
# a.PDB2PDBwLeap()
# a.PDB2FF()
# a.PDBMin()
# a.rm_wat()
# a.PDB2FF()
# a.PDBMD()
