import os
import re
from random import randint
from AmberMaps import *

__doc__='''
This module defines the potential object to operate. Included some common I/O and methods.
-------------------------------------------------------------------------------------
Class PDB
-------------------------------------------------------------------------------------
__init__(self,path)
    Input:  PDB path (any format), MutaFlags
    Output: self.path 
            self.name 
            self.ifformat (Judge by the first line)
          # Connect with the database // add a flow for database?
-------------------------------------------------------------------------------------
Information collecting methods:
-------------------------------------------------------------------------------------
get_tot_resi(self)
    Input:  self (standard Amber format)
    Output: self.tot_resi
get_index_resi(self, index)
get_coord(self, index)
get_seq(self)
get_missing(self, seq)

-------------------------------------------------------------------------------------
PDB operating methods: (changes the self.path to indicated new pdb)
-------------------------------------------------------------------------------------
PDB2PDBwLeap(self, Flag):
PDBMin(self,cycle):
rm_wat(self): remove water and ion for current pdb. (For potential docking)
loopmodel_refine(self): Use different method to model the missing sequence

-------------------------------------------------------------------------------------
Mutation Tools:
-------------------------------------------------------------------------------------
Add_MutaFlag(self,Flag): User assigned Flag or Random Flag using "random"
PDB_check(self):

-------------------------------------------------------------------------------------
Input file generating methods:
-------------------------------------------------------------------------------------
PDB2FF(self):
PDBMD(self):
-------------------------------------------------------------------------------------
'''

class PDB(object):

    path=''
    name=''
    File_str=''
    prmtop_path=''
    inpcrd_path=''
    stage=0 #For debug indicate which stage is the program in
    ifformat=0
    ifcomplete=1

    MutaFlags=[]
    tot_resi=0
    Coord=[]
    raw_sequence={}
    sequence={}

    #current operating residue
    current_index=0
    current_resi_name_3=''
    current_resi_name=''
    current_resi_atom_list=[]
    
    


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
            print('WARNING: This file is not in standard Amber format.')
            print('WARNING: This is normal when refinement/protonation is one of the target steps.\n')


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

    def get_file(self):
        '''
        Get the PDB file and store it in the self.file
        return the self.file
        '''
        self.file = open(self.path)

        return self.file
    
    def get_file_str(self):
        '''
        Get the str of the PDB file and store it in self.File_str
        return the str
        '''
        self.get_file()
        self.File_str =self.file.read()

        return self.File_str



    def get_seq(self):
        '''
        Get sequence for current PDB (self.path) // A general function to obtain the missing residue
        + Re-assign the chain_index base on the order in the file
        + Use "NAN" as a filler to store missing residues (detect internal missing)
        + Check if contain non-standard residue with the Resi_map 
        - Do not include the original residue index
        - Do not include any HETATM (ligand/solvent)
        self.sequence:
        Format: {'Chain_index':['res','res',...]
                 'Chain_index':['res','res',...]
                 ...
                }
        save the info to self.sequence and return it

        self.raw_sequence contain HETATM

        !!NOTE!! the self.sequence needs to be updated after mutation
        '''

        PDB_str = self.get_file_str()

        Chain_str = PDB_str.split('TER')
        
        for chain,i in enumerate(Chain_str):

            Chain_index = chr(65+chain) # Covert to ABC using ACSII mapping
            Chain_sequence=[]
            
            # Get the Chain_sequence
            lines=i.split('\n')

            for line_index, line in enumerate(lines):

                pdb_l = PDB_line(line)

                if pdb_l.line_type == 'ATOM  ' or pdb_l.line_type == 'HETATM':
                    
                    pdb_l.get_resi_index()
                    pdb_l.get_resi_name()
                    
                    # Deal with the first residue
                    if len(Chain_sequence) == 0:
                        Chain_sequence.append(pdb_l.resi_name)
                        last_resi_index = pdb_l.resi_index
                        continue

                    # find next new residue
                    if pdb_l.resi_index != last_resi_index:

                        # Deal with missing residue, fill with "NAN"
                        missing_length = pdb_l.resi_index - last_resi_index - 1 
                        if missing_length > 0:
                            Chain_sequence = Chain_sequence + ['NAN',] * missing_length 
                        
                        # Store the new resi
                        Chain_sequence.append(pdb_l.resi_name)

                    # Update for next loop                
                    last_resi_index = pdb_l.resi_index
            
            self.raw_sequence[Chain_index] = Chain_sequence
        
            self.strip_raw_seq() # strip the raw_sequence and save to sequence

        return self.sequence

    def strip_raw_seq(self):
        '''
        (Used internally) strip the raw_sequence.
        - Delete ligand and solvent
        - Delete chains without residue

        save changes to self.sequence
        '''

        if len(self.raw_sequence) == 0:
            print("The self.raw_sequence should be obtained first")
            raise IndexError

        for chain in self.raw_sequence:
            for resi in self.raw_sequence[chain]:








            










    def PDB2PDBwLeap(self):
        '''
        (The index is correspond to the new index starting from 1 for each chain.)
        Take the 'XA111Y' format for polymer. Re-local the mutation position based on the index of the first residue in the chain (mark: update)

        Use every Flag in self.MutaFlags and tLeap to build mutated structure. 
        **WARNING** if there are multiple mutations on the same index, only the first one will be used.
        Save changes to self.path and self.name
        self.stage=1
        '''

        #Judge if theres same MutaIndex
        for i in range(len(self.MutaFlags)):
            for j in range(len(self.MutaFlags)):
                if i >= j:
                    pass
                else:
                    if self.MutaFlags[i][1] == self.MutaFlags[j][1]:
                        print("**WARNING**")
                        print("There are multiple mutations in the same index, only the first one will be used: "+self.MutaFlags[i][0]+self.MutaFlags[i][1]+self.MutaFlags[i][2])

        # Prepare a label for the filename
        tot_Flag_name=''
        for Flag in self.MutaFlags:
            Flag_name=Flag[0]+Flag[1]+Flag[2]
            tot_Flag_name=tot_Flag_name+'_'+Flag_name

        # Operate the PDB
        out_PDB_path1=self.path[:-4]+tot_Flag_name+'_tmp.pdb'
        out_PDB_path2=self.path[:-4]+tot_Flag_name+'.pdb'

        with open(self.path,'r') as f:
            with open(out_PDB_path1,'w') as of:
                line_index=1
                for line in f:

                    try:
                        match=0
                        
                        for Flag in self.MutaFlags:
                            # Test for every Flag for every lines
                            resi_Index=Flag[1]

                            if line.split()[4] == resi_Index:
                                #This matching method work with ploymer because the standard amber format merge different chains into a same index map
                                #keep if match the OldAtom list

                                match=1
                                
                                #Initialize mutation when match
                                Muta_resi=Flag[2]
                                OldAtoms=['N','H','CA','HA','CB','C','O']
                                #fix for mutate to Gly & Pro
                                if Muta_resi == 'G':
                                    OldAtoms=['N','H','CA','C','O']
                                if Muta_resi == 'P':
                                    OldAtoms=['N','CA','HA','CB','C','O']


                                for i in OldAtoms:
                                    if i == line.split()[2]:
                                        #change the traget residue name !!WARNING!! strong format limitation here!
                                        
                                        new_line=line[:17]+Resi_map[Muta_resi]+line[20:]
                                        
                                        of.write(new_line)
                                        break
                                #Do not keep if not match & in the target residue.       
                                
                                #Dont run for other Flags after first Flag matches. (code: SameSite)
                                break

                        if not match:               
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
        os.system('mv *leap.* *tmp.pdb tleap_cache')

        #Update the file
        self.path = out_PDB_path2
        self.name = self.path.split(os.sep)[-1][:-4]
        self.stage = 1

        return out_PDB_path2


    def Add_MutaFlag(self,Flag):
        '''
        Input: 
        Flags or "random"
        ----------------------------------------------------
        Flag    (e.g. A11B) Can be a str or a list of str.
        ----------------------------------------------------
        Append self.MutaFlags with the Flag.
        (The index is correspond to the new merged index from Standard Amber format for Polymer.)
        ----------------------------------------------------
        'random' or 'r' (default)
        ----------------------------------------------------
        Use self.path & self.sequence // only when ifcomplete
        Save changes appending to self.MutaFlags 
        
        Need to add the chain_index info (mark: update)
        '''

        if type(Flag) == str:

            if Flag == 'r' or Flag == 'random':
                resi_1=''
                resi_2=''
                Muta_idx=''

                # This method takes the first chain of a PDB file to count the total number of residues. Will encounter problem when resi_index of the ligand was inserted between the chain. (mark:enzyme)
                # This method use get_seq and self.sequence. Randomize over the chain and use the 'XA123Y' mutaflag format.  (mark: update)
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
                        
                self.MutaFlags.append((resi_1,Muta_idx,resi_2))

            else:
                self.MutaFlags.append((Flag[0],re.search('[0-9]+',Flag).group(),Flag[-1]))

        if type(Flag) == list:
            for i in Flag:
                self.MutaFlags.append((i[0],re.search('[0-9]+',i).group(),i[-1]))

        print('Current MutaFlags: ',self.MutaFlags)


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
        Flag_name=''
        for Flag in self.MutaFlags:
            sep_Flag_name=Flag[0]+Flag[1]+Flag[2]
            Flag_name=Flag_name+'_'+sep_Flag_name

        out_path=conf_path

        #run sander
        os.system('mpirun -np 8 $AMBERHOME/bin/sander.MPI -O -i '+conf_path+'/min'+tag+'.in -o ' +out_path+'/min' +Flag_name+'.out -p '+self.prmtop_path+' -c '+self.inpcrd_path+' -r '+out_path+'/min'+Flag_name+'.rst -ref '+self.inpcrd_path)
        os.system('mpirun -np 8 $AMBERHOME/bin/sander.MPI -O -i '+conf_path+'/heat'+tag+'.in -o '+out_path+'/heat'+Flag_name+'.out -p '+self.prmtop_path+' -c '+out_path+'/min'+Flag_name+'.rst -ref ' +out_path+'/min'+Flag_name+'.rst -r ' +out_path+'/heat'+Flag_name+'.rst')
        os.system('mpirun -np 8 $AMBERHOME/bin/sander.MPI -O -i '+conf_path+'/equi'+tag+'.in -o '+out_path+'/equi'+Flag_name+'.out -p '+self.prmtop_path+' -c '+out_path+'/heat'+Flag_name+'.rst -ref '+out_path+'/heat'+Flag_name+'.rst -r '+out_path+'/equi'+Flag_name+'.rst')
        os.system('mpirun -np 8 $AMBERHOME/bin/sander.MPI -O -i '+conf_path+'/prod'+tag+'.in -o '+out_path+'/prod'+Flag_name+'.out -p '+self.prmtop_path+' -c '+out_path+'/equi'+Flag_name+'.rst -ref '+out_path+'/equi'+Flag_name+'.rst -r '+out_path+'/prod'+Flag_name+'.rst -x '+out_path+'/prod'+Flag_name+'.nc')

        #Or add a class constant?
        return 'MD/prod'+Flag_name+'.rst'





class PDB_line(object):
    '''
    Class for decoding the line in the PDB file. Functions used internally. 
    '''

    line=''
    line_type=''

    resi_name=''
    resi_index=''


    def __init__(self,line):
        '''
        initilize with a specific line in the PDB file.
        Get self.line_type
        '''
        self.line=line
        self.line_type = self.line[0:6]

    '''
    =====
    Residue
    =====
    '''

    def get_resi_name(self):
        '''
        Get the residue name from self.line
        save the result in self.resi_name
        '''
        self.resi_name = self.line[17:20]

        return self.resi_name
    
    def get_resi_index(self):
        '''
        Get the residue index from self.line
        save the result in self.resi_index
        '''
        self.resi_index = int(self.line[22:26])

        return self.resi_index







# Pre-refinement functions
    def get_missing(self, seq):
        '''
        get_missing(self, seq)
        Compare self.sequence (from the get_seq) with the seq (str from the uniport)
        1. No missing
        2. Terminal missing
        3. Internal missing
        '''
        pass

    def PDB_loopmodel_refine(self, method='Rosetta'):
        '''
        Use different methods to model the missing sequence
        Methods:
        - pyRosetta
        - trRosetta
        '''
        pass

# 在残缺时使用某种占位字母记录？读取时对相邻的序号做差，大于一的添加N-1个空缺符（不能用字符串对比方法重现残缺位置因为确实了链接关系不知道在哪插入）用seq重写随机突变的方法（flag的存储，随机的方式，突变的实施）
# 结合class PDB_line的思路完成序列的提取, 用chr(65-1+i)完成ABC转化

# func outside of the class
def get_PDB(name):
    '''
    connect to the database
    '''
    pass


def PDB_to_AMBER_PDB(path):
    '''
    Make the file convertable with tleap without error
    - Test result: The header part cause duplication of the residues. Deleting that part may give normal tleap output
    - Test result: For some reason, some ligand will miss if directly covert after simple header cutting
    - WARNINGs
    '''
    pass


#TestOnly
# a=PDB(r'1BRS_refine_amber_A150C_A192C_E73F_min_rmW.pdb')
# a.Add_MutaFlag('E73F')
# a.Add_MutaFlag('random')
# a.PDB2PDBwLeap()
# a.PDB2FF()
# a.PDBMin()
# a.rm_wat()
# a.PDB2FF()
# a.PDBMD()
