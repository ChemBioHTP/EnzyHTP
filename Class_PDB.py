import os
import re
import argparse
from random import randint
from AmberMaps import *
from wrapper import *
from Class_Structure import *
from Class_line import *
try:
    from pdb2pqr.main import main_driver as run_pdb2pqr
    from pdb2pqr.main import build_main_parser as build_pdb2pqr_parser
except ImportError:
    raise ImportError('PDB2PQR not installed.')

try:
    import propka.lib
except ImportError:
    raise ImportError('PropKa not installed.')


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
get_if_ligand = get_seq # alias
get_if_art_resi = get_seq # alias
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

class PDB():

    path=''
    name=''
    File_str=''
    prmtop_path=''
    inpcrd_path=''
    pqr_path=''
    stage=0 #For debug indicate which stage is the program in
    
    ifformat=0
    
    if_complete=None
    if_complete_chain={} 
    if_art_resi=None
    if_ligand=None

    MutaFlags=[]
    tot_resi=0
    Coord=[]
    raw_sequence={}
    sequence={}
    sequence_one={}

    #current operating residue
    current_index=0
    current_resi_name_3=''
    current_resi_name=''
    current_resi_atom_list=[]
    
    


    #regex pattern of standard Amber format line
    ATOM_pattern=r'ATOM  [ ,0-9][ ,0-9][ ,0-9][ ,0-9][0-9]  [A-Z][ ,0-9,A-Z][ ,0-9,A-Z] [A-Z][A-Z][A-Z]  [ ,0-9][ ,0-9][ ,0-9][0-9]    [ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9][ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9][ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9](?:         [ ,A-Z][ ,A-Z][A-Z])?'
    HETATM_pattern=r'HETATM[ ,0-9][ ,0-9][ ,0-9][ ,0-9][0-9]  [A-Z][ ,0-9,A-Z][ ,0-9,A-Z] [A-Z][A-Z][A-Z]  [ ,0-9][ ,0-9][ ,0-9][0-9]    [ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9][ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9][ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9](?:         [ ,A-Z][ ,A-Z][A-Z])?'


    def __init__(self,PDB_PATH='',PDB_File_str=''):
        '''
        Waiting for reform
        '''
        #initilize
        self.path=''
        self.name=''
        self.File_str=''
        self.prmtop_path=''
        self.inpcrd_path=''
        self.pqr_path=''
        self.stage=0         
        self.ifformat=0        
        self.if_complete=None
        self.if_complete_chain={} 
        self.if_art_resi=None
        self.if_ligand=None
        self.MutaFlags=[]
        self.tot_resi=0
        self.Coord=[]
        self.raw_sequence={}
        self.sequence={}
        self.sequence_one={}
        self.current_index=0
        self.current_resi_name_3=''
        self.current_resi_name=''
        self.current_resi_atom_list=[]
        #initilize


        if len(PDB_PATH) == 0:
            self.File_str=PDB_File_str
            # 存一个文件和path？
        else:
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
                pass
                #print('WARNING: This file is not in standard Amber format.')
                #print('WARNING: This is normal when refinement/protonation is one of the target steps.\n')


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
        Get the PDB file and store it in the self.File
        return the self.File
        '''
        self.File = open(self.path)

        return self.File
    
    def get_file_str(self):
        '''
        Get the str of the PDB file and store it in self.File_str
        return the str
        '''
        if self.File_str == '':
            self.get_file()
            self.File_str =self.File.read()

        return self.File_str


    '''
    =========
    Sequence (not require Amber format)
    =========
    '''

    def get_seq(self, Oneletter=0):
        '''
        get_seq(self, Oneletter=0)
        Support most PDB types
        ----------------------------
        Oneletter
        = 0 // use 3-letter format to represet each residue
        = 1 // use 1-letter format to represet each residue
        ----------------------------
        Get sequence for current PDB (self.path) // A general function to obtain the missing residue
        + Use "NAN"/"-" as a filler to store missing residues (detect internal missing)
        + Check if contain non-standard residue with the Resi_map --> self.if_art_resi
        + Check if contain ligand with the Resi_map and TIP3P_map --> self.if_ligand
        - (WARNING) Require the ligand in seperate chains
        - Re-assign the chain_index base on the order in the file
        - Do not include the original residue index
        - Do not include any HETATM (ligand/solvent)
        ----------------------------
        save the info to self.sequence/self.sequence_one and return it
        ----------------------------
        self.sequence:
        Format: {'Chain_index':['res','res',...]
                 'Chain_index':['res','res',...]
                 ...
                }

        self.raw_sequence:
            Internal used, containing HETATM

        self.sequence_one: (when Oneletter=1)
        Format: {'Chain_index':'seq'
                 'Chain_index':'seq'
                 ...
                }

        !!NOTE!! the self.sequence needs to be updated after mutation
        '''

        PDB_str = self.get_file_str()

        Chain_str = PDB_str.split('\nTER') # Note LF is required
        
        for chain,i in enumerate(Chain_str):

            Chain_index = chr(65+chain) # Covert to ABC using ACSII mapping
            Chain_sequence=[]
            
            # Get the Chain_sequence
            lines=i.split('\n')

            for line in lines:

                pdb_l = PDB_line(line)

                if pdb_l.line_type == 'ATOM' or pdb_l.line_type == 'HETATM':
                    
                    # Deal with the first residue
                    if len(Chain_sequence) == 0:
                        Chain_sequence.append(pdb_l.resi_name)
                        last_resi_index = pdb_l.resi_id
                        continue

                    # find next new residue
                    if pdb_l.resi_id != last_resi_index:

                        # Deal with missing residue, fill with "NAN"
                        missing_length = pdb_l.resi_id - last_resi_index - 1 
                        if missing_length > 0:
                            Chain_sequence = Chain_sequence + ['NAN',] * missing_length 
                        
                        # Store the new resi
                        Chain_sequence.append(pdb_l.resi_name)

                    # Update for next loop                
                    last_resi_index = pdb_l.resi_id
            
            self.raw_sequence[Chain_index] = Chain_sequence
            
        self._strip_raw_seq() # strip the raw_sequence and save to sequence
        
        self.get_if_complete()

        if Oneletter == 1:
            self._get_Oneletter()
            return self.sequence_one
        else:
            return self.sequence

    get_if_ligand = get_seq # alias
    get_if_art_resi = get_seq

    def _strip_raw_seq(self):
        '''
        (Used internally) strip the raw_sequence.
        - Delete ligand and solvent
        - Delete chains without residue

        save changes to self.sequence

        Judge if containing any ligand or artificial residue

        save to self.if_ligand and self.if_art_resi
        '''
        new_index=0

        #if the value is not changed to 1 then it's 0
        self.if_art_resi = 0
        self.if_ligand = 0

        if len(self.raw_sequence) == 0:
            print("The self.raw_sequence should be obtained first")
            raise IndexError

        for chain in self.raw_sequence:
            
            chain_seq=[]
            if_realchain = 0

            # Judge if real chain
            for name in self.raw_sequence[chain]:
                clean_name = name.strip(' ')
                if clean_name in Resi_map2.keys():
                    if_realchain = 1
            
            if if_realchain:
                for name in self.raw_sequence[chain]:
                    clean_name = name.strip(' ')                               
                    chain_seq.append(clean_name)

                    # An artificial residue will be a residue in a realchain but not included in force field map
                    if clean_name not in Resi_map2 and clean_name != 'NAN':
                        self.if_art_resi = 1                        
                        ## PLACE HOLDER for further operation on artificial residue ##

            else:
                # Judge if containing any ligand
                for name in self.raw_sequence[chain]:
                    clean_name = name.strip(' ') 
                    if clean_name not in TIP3P_map and clean_name != 'NAN':
                        self.if_ligand = 1
                        ## PLACE HOLDER for further operation on artificial residue ##
            
            # only add realchain to the self.sequence
            if len(chain_seq) != 0:
                chain_Index=chr(new_index+65)
                self.sequence[chain_Index]=chain_seq
                new_index=new_index+1


    def _get_Oneletter(self):
        '''
        (Used internally) convert sequences in self.sequence to oneletter-based str
        - The 'NAN' is convert to '-'
        - Present unnature residue as full 3-letter name

        save to self.sequence_one
        '''
        if len(self.sequence) == 0:
            print("The self.sequence should be obtained first")
            raise IndexError

        for chain in self.sequence:
            chain_Seq=''
            for name in self.sequence[chain]:
                c_name=name.strip(' ')
                if c_name == 'NAN':
                    chain_Seq=chain_Seq+'-'
                else:
                    if c_name in Resi_map2:
                        chain_Seq=chain_Seq+Resi_map2[c_name]
                    else:
                        chain_Seq=chain_Seq+' '+c_name+' '
                
            self.sequence_one[chain]=chain_Seq


    def get_if_complete(self):
        '''
        Judge if the self.sequence (from the get_seq) has internal missing parts
        Save the result to:
            self.if_complete ------- for intire PDB
            self.if_complete_chain - for each chain
        '''
        if len(self.sequence.keys()) == 0:
            print('Please get the sequence first')
            raise IndexError
        
        self.if_complete=1 # not flow in if then 1
        
        for chain in self.sequence:
            self.if_complete_chain[chain]=1 # not flow in if then 1
            for resi in self.sequence[chain]:
                if resi == 'NAN':
                    self.if_complete=0
                    self.if_complete_chain[chain]=0
                    break
        return self.if_complete, self.if_complete_chain


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

    '''
    ========
    Protonation
    ========
    '''
    def get_protonation(self, ph=7.0):
        '''
        Get protonation state based on PDB2PQR:
        1. Use PDB2PQR, save output to self.pqr_path
        2. Fix problems:
            - Metal center: 
            (detect donor(base on atom type) in a metal type based radiis) --> open for customize
                - Fix1: deprotonate all. 
                - Fix2: rotate if there're still lone pair left 
                - Fix3: run pka calculate containing ion (maybe pypka) and run Fix2 based on the result
                # switch HIE HID when dealing with HIS
        save to self.path
        '''
        out_path=self.name+'_aH.pdb'
        self._get_protonation_pdb2pqr(ph=ph)
        self._protonation_Fix(out_path)

 
    def _get_protonation_pdb2pqr(self,ffout='AMBER',ph=7.0,out_path=''):
        '''
        Use PDB2PQR to get the protonation state for current PDB. (self.path)
        current implementation just use the outer layer of PDB2PQR. Update to inner one and get more infomation in the furture.
        
        save the result to self.pqr_path
        '''
        # set default value for output pqr path
        if len(out_path) == 0:
            self.pqr_path = self.name+'.pqr'
        else:
            self.pqr_path = out_path
        
        # input of PDB2PQR
        pdb2pqr_parser = build_pdb2pqr_parser()
        args = pdb2pqr_parser.parse_args(['--ff=PARSE','--ffout='+ffout,'--with-ph='+str(ph),self.path,self.pqr_path])
        # use context manager to hide output
        with HiddenPrints('./._get_protonation_pdb2pqr.log'):
            run_pdb2pqr(args)


    def _protonation_Fix(self, out_path, Metal_Fix='1'):
        '''
        Add in the missing atoms and run detailed fixing
        save to self.path
        '''

        # Add missing atom (from the PDB2PQR step. Update to func result after update the _get_protonation_pdb2pqr func)       
        # Now metal only
        stru = structure(self)
        # find Metal center and combine with the pqr file
        metal_list = stru.metalatoms
        with open(self.pqr_path) as f:
            lines=f.readlines()
            new_lines=[]
            for metal in metal_list:
                new_lines.append(metal.get_line(ff='AMBER'))
            new_lines.append('TER\n')
            lines = lines[:-2] + new_lines + lines[-1]
            with open(out_path,'w') as of:
                for line in lines:
                    of.write(line)
        self.path = out_path

        # fix metal donor
        self._protonation_Metal_Fix(Fix=Metal_Fix)        

    def _protonation_Metal_Fix(self, Fix):
        '''
        fix the protonation state of donor residue of the metal center.
        overwrite the self.path
        '''
        stru = structure(self)
        metal_list=stru.find_metal()

        if Fix == '1':
            for metal in metal_list:
                d_resi_list = metal.get_donor_residue()
                for d_resi in d_resi_list:
                    d_resi.deprotonate()
            stru.build(ff='AMBER')

        if Fix == '2':
            pass

        if Fix == '3':
            pass


    '''
    ========
    Mutation
    ========
    '''


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
        Use self.path & self.sequence // only when if_complete
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


    '''
    ========
    Gerneral MD
    ========
    '''

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
        '''
        Run a minization use self.prmtop and self.inpcrd
        Save changed PDB to self.path (containing water and ions)
        '''
        
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
        Now only skip [Na+,Cl-,WAT,HOH] // Append more in the future.
        Save changed files into self.path.
        '''
        out_path = self.path[:-4]+'_rmW.pdb'

        with open(self.path) as f:
            with open(out_path,'w') as of:
                skip_list=['Na+','Cl-','WAT','HOH']
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









# 主函数的metal处理已经写好
# 之后再structure类对应主函数的使用写具体的函数

# 1B8Q这个文件有问题，待处理此类
# 了解并解决加氢的问题
# 现在用pdb2pqr解决了一部分问题：
# - 改了其中的两个bug
# - 但是金属中心相关的问题还没有完成 (写一个单独的函数删去加多余的H) # htmd保留了离子但没有很好的处理周围的残基的质子态 # pypka貌似可以带离子计算

#pypka
#貌似必须要amber处理之前的结构，如果变成HIE等名称将无法识别

# 判断库里有多少可以直接跑, 遗留问题：无法区分共结晶用的一些小分子和配体。处理酶的时候一概删去
# 写pdb转化函数：
# - 配体存储于独立的链 (貌似都是，也可以判断一下，对于所有判断为人工残基的)
# - 可以被leap正常读取

# 跑这部分的MD，重写conf类
# 
# 基于seq完成get_missing的判断，
# 写修复的方法Rosetta优先 // 与uniport的对比
#
# 用seq重写随机突变的方法（flag的存储，随机的方式，突变的实施）
#
# 尝试适配直接传入字符串？或许可以写一个保存文件的函数? （待定）还是要以外部文件为主，因为要和别的程序交互

#htmd的不足
#1.结构修复 （只做到保留残基序号从而保留缺失片段的位置）
#2.金属中心附近的质子态 （可以保留金属离子，但是质子态处理的很差） 

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
    - Test result: `cat 1NVG.pdb |grep "^ATOM\|^HETATM\|^TER|^END" > 1NVG-grep.pdb` will be fine
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
