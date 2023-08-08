import copy
from math import ceil
import os
import io
import re
import pandas as pd
import shutil
from shutil import rmtree
from subprocess import SubprocessError, run, CalledProcessError
from random import choice
from typing import Dict, Union, List
from AmberMaps import *
from wrapper import *
from Class_Structure import *
from Class_line import *
from Class_Conf import Config, Layer
from Class_ONIOM_Frame import *
from core import job_manager
from core.clusters._interface import ClusterInterface
from helper import (
    Conformer_Gen_wRDKit, 
    delete_idx_line, 
    decode_atom_mask, 
    get_center, 
    get_field_strength_value, 
    line_feed, 
    mkdir, 
    generate_Rosetta_params,
    AmberError,
    )
try:
    from pdb2pqr.main import main_driver as run_pdb2pqr
    from pdb2pqr.main import build_main_parser as build_pdb2pqr_parser
except ImportError:
    raise ImportError('PDB2PQR not installed.')

try:
    import propka.lib
except ImportError:
    raise ImportError('PropKa not installed.')
try:
    import openbabel
    import openbabel.pybel as pybel
except ImportError:
    raise ImportError('OpenBabel not installed.')



__doc__='''
This module handles file I/O and external programs.
-------------------------------------------------------------------------------------
Class PDB
-------------------------------------------------------------------------------------
__init__(self, PDB_input, wk_dir = '', name = '', input_type='path')
-------------------------------------------------------------------------------------
Information collecting methods:
-------------------------------------------------------------------------------------
get_stru(self, ligand_list=None)
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
>>>>>>> Develop Note <<<<<<<<<
For every potential first use of self.file_str or self.path or self.stru
Will do nothing if the correct var is already exist.
 - use _get_file_str
 - use _get_file_path
 - use get_stru 
'''

class PDB():
    def __init__(self, PDB_input, wk_dir = '', name = '', input_type='path'):
        '''
        initiate PDB object
        -------------------
        PDB_input
        wk_dir      : working directory (default: current dir ) **all file in the workflow are constructed base on this.**
        name        : self assigned filename (default: from path or UNKNOW if no path)
        input_type  : path (default) / file_str / file
        -------------------
        '''
        # Nessessary initilize for empty judging
        self.stru = None
        self.path = None
        self.prmtop_path = None
        self.MutaFlags = []
        self.nc=None
        self.frames=None
        self.prepi_path = {}
        self.frcmod_path = {}
        self.disulfied_residue_pairs = []
        # default MD conf.
        self._init_MD_conf()
        # default ONIOM layer setting
        self.layer_preset = Config.Gaussian.layer_preset
        self.layer_atoms = Config.Gaussian.layer_atoms

        if wk_dir == '':
            #current dir by default
            self.dir = '.'
        else:
            if wk_dir[-1] == '/':
                self.dir = wk_dir[:-1]
            else:
                self.dir = wk_dir

        if input_type == 'path':
            self.path = PDB_input
            self._update_name()
        if input_type == 'file_str':
            self.file_str = PDB_input
            if name == '':
                # default name
                name = 'UNKNOW'
            self.name = name
            self.path_name = self.dir+'/'+self.name
        if input_type == 'file':
            self.file_str = PDB_input.read()
            if name == '':
                # default name
                name = 'UNKNOW'
            self.name = name
            self.path_name = self.dir+'/'+self.name
        if input_type not in ['path', 'file_str', 'file']:
            raise Exception('PDB.__init__: only take path or file_str or file.')

        # make cache
        self.cache_path = self.dir+'/cache'
        mkdir(self.cache_path)
        

    def _update_name(self):
        '''
        update name
        '''
        suffix_len = len(self.path.split('.')[-1]) + 1
        self.name=self.path.split(os.sep)[-1][:-suffix_len]
        self.path_name = self.dir+'/'+self.name
        

    def get_stru(self, ligand_list=None, renew = 0):
        '''
        Convert current PDB file (self.path) to a Struture object (self.stru).
        ------
        input_name  : a name tag for the object (self.name by default)
        ligand_list : a list of user assigned ligand residue names.
        renew       : 1: force generate a new stru_obj
        '''
        # indicated by self.name
        input_name = self.name
        # if get new stru
        get_flag = 0
        if self.stru is not None:
            if self.stru.name != self.name:
                get_flag = 1
                # warn if possible wrong self.stru
                if Config.debug >= 1:
                    print('PDB.get_stru: WARNING: self.stru has a different name')
                    print('     -self.name: '+self.name)
                    print('     -self.stru.name: '+self.stru.name)
        else:
            get_flag = 1

        if get_flag or renew:
            if Config.debug >= 1:
                print('PDB.get_stru: Getting new stru')
            if self.path is not None:
                self.stru = Structure.fromPDB(self.path, input_name=input_name, ligand_list=ligand_list)
            else:
                self.stru = Structure.fromPDB(self.file_str, input_type='file_str', input_name=input_name, ligand_list=ligand_list)
        elif Config.debug >= 1:
            print('PDB.get_stru: Not getting new stru - have existing self.stru with the same name and renew == 0')

    def _get_file_str(self):
        '''
        read file_str from path is not pre-exist.
        -------------
        recommend before every potential first use of self.file_str
        '''
        if self.path is None:
            return self.file_str
        else:
            return open(self.path).read()


    def _get_file_path(self):
        '''
        save a file and get path if self.path is None
        -------------
        recommend before every potential first use of self.path
        '''
        if self.path is None:
            self.path = self.path_name+'.pdb'
            with open(self.path,'w') as of:
                of.write(self.file_str)
        return self.path


    def _init_MD_conf(self):
        '''
        initialize default MD configuration. Replaced by manual setting if assigned later. TODO may be need copy & operation in where it is used.
        '''
        self.conf_min = Config.Amber.conf_min
        self.conf_heat = Config.Amber.conf_heat
        self.conf_equi = Config.Amber.conf_equi
        self.conf_prod = Config.Amber.conf_prod


    def set_oniom_layer(self, atom_list=[], preset=0):
        '''
        set oniom layer with san check
        layer_atoms are in higher pirority
        '''
        if len(atom_list) != 0:
            if all(type(x) is str for x in atom_list):
                self.layer_atoms = atom_list
            else:
                raise Exception('set_oniom_layer: atom_list require a list of str. e.g.: [\'1-9\',\'11,13-15\']')
        else:
            if preset == 0:
                raise Exception('set_oniom_layer: please assign one of the argument: atom_list (a list of layer_atoms) or preset (other than 0)')
            else:
                self.layer_preset = preset


    def get_last_A_id(self):
        '''
        get last atom id in PDB
        '''
        with open(self.path) as f:
            fl = f.readlines()
            for i in range(len(fl)-1,-1,-1):
                pdbl = PDB_line(fl[i])
                if pdbl.line_type == 'ATOM' or pdbl.line_type == 'HETATM':
                    return pdbl.atom_id
    '''
    =========
    Sequence TODO: reform(most of it has been moved to class Chain，only judgement of art residue and overlook of whole structure remains)
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

        PDB_str = self._get_file_str()
        self.raw_sequence = {}
        self.sequence = {}
        self.sequence_one = {}

        Chain_str = PDB_str.split(line_feed+'TER') # Note LF is required
        
        for chain,i in enumerate(Chain_str):

            Chain_index = chr(65+chain) # Covert to ABC using ACSII mapping
            Chain_sequence=[]
            
            # Get the Chain_sequence
            lines=i.split(line_feed)

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
        
        self.if_complete=1 # if not flow in then 1
        self.if_complete_chain = {}
        
        for chain in self.sequence:
            self.if_complete_chain[chain]=1 # if not flow in then 1
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
        remember to update the name.
        '''
        pass

    '''
    ========
    Protonation
    ========
    '''
    def get_protonation(self, ph=7.0, keep_id=0, if_prt_ligand=1):
        '''
        Get protonation state based on PDB2PQR:
        1. Use PDB2PQR, save output to self.pqr_path
        2. Fix problems:
            - Metal center: 
            (detect donor(base on atom type) in a metal type based radiis) --> open for customize
                - Fix1: deprotonate all. 
                - Fix2: rotate if there're still lone pair left 
                - Fix3: run pka calculate containing ion (maybe pypka) and run Fix2 based on the result
            - Ligand:
                - Use OpenBable to protonate ligand by default
                # switch HIE HID when dealing with HIS
        save to self.path
        ----------
        ph: pH when determine the protonation state
        keep_id: if keep ids of original pdb file
        if_prt_ligand: if re-protonate ligand. (since sometime its already protonated)
        '''
        out_path=self.path_name+'_aH.pdb'
        self._get_file_path()
        disulfied_residue_pairs = self._get_protonation_pdb2pqr(ph=ph)
        self.disulfied_residue_pairs = disulfied_residue_pairs # TODO this is a temp solution. This requires runnning this function everytime after mutation which is also recommanded.
        self._protonation_Fix(out_path, ph=ph, keep_id=keep_id, if_prt_ligand=if_prt_ligand)
        self.path = out_path
        self._update_name()
        self.stru.name=self.name


    def _get_protonation_pdb2pqr(self,ffout='AMBER',ph=7.0,out_path=''):
        '''
        Use PDB2PQR to get the protonation state for current PDB. (self.path)
        current implementation just use the outer layer of PDB2PQR. Update to inner one and get more infomation in the furture. 
            (TARGET: 1. what is deleted from the structure // metal, ligand)
        
        save the result to self.pqr_path

        Return:
            return disulfied_residue_pairs in a format of e.g.: [(('A', 496), ('A', 440)), ...]
        '''
        # set default value for output pqr path
        if len(out_path) == 0:
            self.pqr_path = self.path_name+'.pqr'
        else:
            self.pqr_path = out_path
        
        # input of PDB2PQR
        pdb2pqr_parser = build_pdb2pqr_parser()
        pdb2pqr_logging_level = 'INFO'
        args = pdb2pqr_parser.parse_args(['--ff=PARSE',f'--log-level={pdb2pqr_logging_level}','--ffout='+ffout,'--with-ph='+str(ph),self.path,self.pqr_path])
        # use context manager to hide output
        with HiddenPrints(f'{self.dir}/._get_protonation_pdb2pqr.log'):
            missed_residues, pka_df, biomolecule = run_pdb2pqr(args)

        # S-S bond
        disulfied_residue_pairs = []
        for res in biomolecule.residues:
            if getattr(res, 'ss_bonded', None):
                res_key = (res.chain_id, res.res_seq)
                for exist_pair in disulfied_residue_pairs:
                    if res_key in exist_pair:
                        break
                else:
                    disulfied_residue_pairs.append(
                        (res_key,
                        (res.ss_bonded_partner.chain_id, res.ss_bonded_partner.res_seq)))
        for ss_bond in disulfied_residue_pairs:
            print(f'INFO: detected disulfied bond by PDB2PQR: {ss_bond[0]} - {ss_bond[1]}')

        return disulfied_residue_pairs


    def _protonation_Fix(self, out_path, Metal_Fix='1', ph = 7.0, keep_id=0, if_prt_ligand=1):
        '''
        Add in the missing atoms and run detailed fixing
        save to self.path
        '''

        # Add missing atom (from the PDB2PQR step. Update to func result after update the _get_protonation_pdb2pqr func)       
        # Now metal and ligand

        old_stru = Structure.fromPDB(self.path)
        new_stru = Structure.fromPDB(self.pqr_path)

        # find Metal center and combine with the pqr file
        metal_list = old_stru.get_metal_center()
        if len(metal_list) > 0:
            new_stru.add(metal_list, sort = 0)
            # fix metal environment
            new_stru.protonation_metal_fix(Fix = 1)

        # protonate ligands and combine with the pqr file
        if len(old_stru.ligands) > 0:
            lig_dir = self.dir+'/ligands/'
            mkdir(lig_dir)
            lig_paths = [(i,k) for i,j,k in old_stru.build_ligands(lig_dir, ifname=1)]

            new_ligs = []
            for lig_path, lig_name in lig_paths:
                if if_prt_ligand:
                    new_lig_path, net_charge = self.protonate_ligand(lig_path, ph=ph)
                else:
                    # keep original structure
                    new_lig_path, net_charge = (lig_path, None)
                new_ligs.append(Ligand.fromPDB(new_lig_path, resi_name=lig_name, net_charge=net_charge, input_type='path'))
            new_stru.add(new_ligs, sort = 0)

        # PLACE HOLDER for other fix

        # build file
        if not keep_id:
            new_stru.sort()
        new_stru.build(out_path, keep_id=keep_id)
        self.stru = new_stru


    @classmethod
    def protonate_ligand(cls, path, method='PYBEL', ph = 7.0, keep_name=1):
        '''
        Protonate the ligand from 'path' with 'method', provide out_path and net charge.
        TODO "obabel -ipdb ligand_1.pdb -opdb pdb -O ligand_1_aHt.pdb -h" can keep names, but how is it accessed by pybel
        ---------------
        method      : PYBEL (default)
                      Dimorphite (from https://durrantlab.pitt.edu/dimorphite-dl/) TODO seems better and with better python API.
                      OPENBABEL (not working if block warning output)
        ph          : 7.0 by default 
        keep_name   : if keep original atom names of ligands (default: 1)
                        - check if there're duplicated names, add suffix if are.
        '''
        outp1_path = path[:-4]+'_badname_aH.pdb'
        out_path = path[:-4]+'_aH.pdb'
        # outm2_path = path[:-4]+'_aH.mol2'

        if method == 'OPENBABEL':
            # not working if block warning output for some reason
            # openbabel.obErrorLog.SetOutputLevel(0)
            obConversion = openbabel.OBConversion()
            obConversion.SetInAndOutFormats("pdb", "pdb")
            mol = openbabel.OBMol()
            obConversion.ReadFile(mol, path)
            mol.AddHydrogens(False, True, ph)
            obConversion.WriteFile(mol, out_path)
        if method == 'PYBEL':
            pybel.ob.obErrorLog.SetOutputLevel(0)
            mol = next(pybel.readfile('pdb', path))
            mol.OBMol.AddHydrogens(False, True, ph)
            mol.write('pdb', outp1_path, overwrite=True)
            # fix atom label abd determing net charge
            if keep_name:
                cls._fix_ob_output(outp1_path, out_path, ref_name_path=path)
            else:
                cls._fix_ob_output(outp1_path, out_path)
            # determine partial charge
            # > METHOD 1<
            net_charge = cls._ob_pdb_charge(outp1_path)
            # > METHOD 2 <
            # mol.write('mol2', outm2_path, overwrite=True)
            # mol = next(pybel.readfile('mol2', outm2_path))
            # net_charge=0
            # for atom in mol:
            #     net_charge=net_charge+atom.formalcharge
        if method == 'Dimorphite':
            pass
        return out_path, net_charge


    @classmethod
    def _fix_ob_output(cls, pdb_path, out_path, ref_name_path=None):
        '''
        fix atom label in pdb_pat write to out_path
        ---------
        ref_name_path: if use original atom names from pdb
        - default: None
            according to tleap output, the name could be just *counting* the element start from ' ' to number
        - : not None
            check if there're duplicated names originally, add suffix if there are.
        '''
        if ref_name_path != None:
            ref_a_names = []
            with open(ref_name_path) as rf:
                pdb_ls = PDB_line.fromlines(rf.read())
                ref_resi_name = pdb_ls[0].resi_name
                for pdb_l in pdb_ls:
                    if pdb_l.line_type == 'HETATM' or pdb_l.line_type == 'ATOM':
                        # pybel use line order (not atom id) to assign new atom id
                        ref_a_names.append(pdb_l.atom_name)

        with open(pdb_path) as f:
            with open(out_path, 'w') as of:
                # count element in a dict
                ele_count={}
                pdb_ls = PDB_line.fromlines(f.read())
                line_count = 0
                for pdb_l in pdb_ls:
                    if pdb_l.line_type == 'HETATM' or pdb_l.line_type == 'ATOM':
                        if ref_name_path == None:
                            ele = pdb_l.get_element()
                        else:
                            if line_count < len(ref_a_names):    
                                ele = ref_a_names[line_count]
                            else:
                                ele = pdb_l.get_element() # New atoms
                            pdb_l.resi_name = ref_resi_name
                            line_count += 1
                        # determine the element count
                        try:
                            # rename if more than one (add count)
                            ele_count[ele] += 1
                            pdb_l.atom_name = ele+str(ele_count[ele])
                        except KeyError:
                            ele_count[ele] = 0
                            pdb_l.atom_name = ele
                        of.write(pdb_l.build())
            

    @classmethod
    def _ob_pdb_charge(cls, pdb_path):
        '''
        extract net charge from openbabel exported pdb file
        '''
        with open(pdb_path) as f:
            net_charge=0
            pdb_ls = PDB_line.fromlines(f.read())
            for pdb_l in pdb_ls:
                if pdb_l.line_type == 'HETATM' or pdb_l.line_type == 'ATOM':
                    if len(pdb_l.get_charge()) != 0:
                        charge = pdb_l.charge[::-1]
                        if Config.debug > 1:
                            print('Found formal charge: '+pdb_l.atom_name+' '+charge)
                        net_charge = net_charge + int(charge)
            return net_charge


    '''
    ========
    Docking
    ========
    '''
    def dock_reactive_substrate(self, substrate, reactive_define, local_lig = 0):
        '''
        Dock substrates into the apo enzyme in a reactive conformation
        Update self.path after docking. Will not update self.path in the middle of the docking
        -----
        In:     "apo-enzyme structure" self.path (pdb)
                "substrate structure" substrate (smile/sdf)
                "constraint parameters" reactive_define - a set of constraint information that defines *Reactive*
                    [reacting residue]
                        id
                        atom_name_1
                    [substrate]
                        atom_name_1
                    --- Will Generate ---
                    [reacting residue]
                        name 
                        atom_name_2
                        atom_name_3
                    [substrate]
                        id 
                        name 
                        atom_name_2 
                        atom_name_3
                    [distance]
                        ideal
                        allowed deviation
                        force constant
                        if_covalent
                    [angle]
                        ideal
                        allowed deviation
                        force constant
        Out:    Enzyme-Substrate Reactive Complex
        -----
        '''
        apo_enzyme = self.path
        if local_lig:
            lig_dir = self.dir+'/ligands/'
            met_dir = self.dir+'/metalcenters/'
        else:
            lig_dir = self.dir+'/../ligands/'
            met_dir = self.dir+'/../metalcenters/'
        mkdir(lig_dir)
        mkdir(met_dir)

        # Build ligands 
        cofactors = self.stru.build_ligands(lig_dir, ifcharge=1, ifunique=1)
        # cofactor_params, cofactors_pdb = generate_Rosetta_params(cofactors, lig_dir, resn='same', out_pdb_name='same')
        # # clean enzyme
        # self.fix_residue_names_for_rosetta()
        # # prepare substrate
        # sub_confs = Conformer_Gen(substrate, method='rdkit')
        # sub_params, sub_pdb, sub_confs_pdb = generate_Rosetta_params(sub_confs, lig_dir, resn='SUB', out_pdb_name='SUB')
        # # define reactive
        # -- search for RULE OF INPUT ID and get a map --
        # self.RDock_generate_cst_file(reactive_define, name_map?)        
        # # initial placement
        # self._RDock_initial_placement(sub_pdb)
        # self._RDock_add_remark_line()
        # # make config file for Rosetta
        # score_path = XXX
        # out_pdb_lib_path = XXX # control the output
        # option_path = Config.Rosetta.generate_xmlNoption(self.path, [cofactor_params, sub_params], score_path, nstruct=100, task='RDock')
        # # Run rosetta script
        # PDB.Run_Rosetta_Script(option_path)
        # Final_pdbs = self.RDock_get_result(score_path, out_pdb_lib_path)
        # self.Rosetta_pdb_to_Amber()



    '''
    ========
    Mutation
    ========
    '''

    def PDB2PDBwLeap(self):
        '''
        Apply mutations using tleap. Save mutated structure PDB in self.path
        ------------------------------
        Use MutaFlag in self.MutaFlags
        Grammer (from Add_MutaFlag):
        X : Original residue name. Leave X if unknow. 
            Only used for build filenames. **Do not affect any calculation.**
        A : Chain index. Determine by 'TER' marks in the PDB file. (Do not consider chain_indexs in the original file.)
        11: Residue index. Strictly correponding residue indexes in the original file. (NO sort applied)
        Y : Target residue name.
        * will just load and save in tleap if MutaFlag is only 'WT'.

        **WARNING** if there are multiple mutations on the same index, only the first one will be used.
        '''

        #Judge if there are same MutaIndex (c_id + r_id)
        for i in range(len(self.MutaFlags)):
            for j in range(len(self.MutaFlags)):
                if i >= j:
                    pass
                else:
                    if (self.MutaFlags[i][1], self.MutaFlags[i][2]) == (self.MutaFlags[j][1], self.MutaFlags[j][2]):
                        if Config.debug >= 1:
                            print("PDB2PDBwLeap: There are multiple mutations at the same index, only the first one will be used: "+self.MutaFlags[i][0]+self.MutaFlags[i][1]+self.MutaFlags[i][2])

        # Prepare a label for the filename
        tot_Flag_name=''
        for Flag in self.MutaFlags:
            Flag_name=self._build_MutaName(Flag)
            tot_Flag_name=tot_Flag_name+'_'+Flag_name

        # Operate the PDB
        out_PDB_path1=self.cache_path+'/'+self.name+tot_Flag_name+'_tmp.pdb'
        out_PDB_path2=self.path_name+tot_Flag_name+'.pdb'

        self._get_file_path()
        with open(self.path,'r') as f:
            with open(out_PDB_path1,'w') as of:
                chain_count = 1
                for line in f:
                    pdb_l = PDB_line(line)
                                        
                    TER_flag = 0
                    if pdb_l.line_type == 'TER':
                        TER_flag = 1
                    # add chain count in next loop for next line
                    if TER_flag:
                        chain_count += 1
                    match=0
                    # only match in the dataline and keep all non data lines
                    if pdb_l.line_type == 'ATOM':
                        for Flag in self.MutaFlags:
                            if 'WT' in Flag:
                                continue
                            # Test for every Flag for every lines
                            t_chain_id=Flag[1]
                            t_resi_id =Flag[2]

                            if chr(64+chain_count) == t_chain_id:
                                if pdb_l.resi_id == int(t_resi_id):
                                    
                                    # do not write old line if match a MutaFlag
                                    match=1
                                    # Keep OldAtoms of targeted old residue
                                    resi_2 = Flag[3]
                                    OldAtoms=['N','H','CA','HA','CB','C','O']
                                    #fix for mutations of Gly & Pro
                                    if resi_2 == 'G':
                                        OldAtoms=['N','H','CA','C','O']
                                    if resi_2 == 'P':
                                        OldAtoms=['N','CA','HA','CB','C','O']

                                    for i in OldAtoms:
                                        if i == pdb_l.atom_name:
                                            new_line=line[:17]+Resi_map[resi_2]+line[20:]                                        
                                            of.write(new_line)
                                            break                                
                                    #Dont run for other Flags after first Flag matches.
                                    break

                    if not match:               
                        of.write(line)


        # Run tLeap 
        #make input
        leapin_path = self.cache_path+'/leap_P2PwL.in'
        leap_input=open(leapin_path,'w')
        leap_input.write('source leaprc.protein.ff14SB\n')
        leap_input.write('a = loadpdb '+out_PDB_path1+'\n')
        leap_input.write('savepdb a '+out_PDB_path2+'\n')
        leap_input.write('quit\n')
        leap_input.close()
        #run
        os.system('tleap -s -f '+leapin_path+' > '+self.cache_path+'/leap_P2PwL.out')
        if Config.debug <= 1:
            os.system('rm leap.log')

        #Update the file
        self.path = out_PDB_path2
        self._update_name()

        return self.path


    def Add_MutaFlag(self, Flag : str = 'r', if_U : bool = 0, if_self : bool = 0):
        """Determine which mutation to deploy to the structure.
        
        User can 1. assign specific mutation(s). 
              or 2. assign random mutation(s) by specifing a rule of randomlization.
        The function will 
        1. decode the input to a list of python tuples that 
            each contains the information of a mutation as ('X','A','###','Y') 
            save to *self.MutaFlags*
        2. check the if the assigned mutation is reasonable that
            specifing the correct original residue
            within the range of chain ids
            within the range of residue ids
            the target residue is within the range of canonical AAs

        Args:
        (self):
        The requirements for the input pdb object are:
            TODO(shaoqz): add after refactoring
        Flag: 
            str or list of strings indicating specifc mutation(s) or a rule of randomlization. (default: r)
            Grammer:
            **Assign specific mutation(s)**
                'XA##Y'
                For each str that represent a mutation, the format should be 'XA##Y', in which
                X : Original residue letter (One-letter code). Leave X if unknow. 
                    (Only used for checking that will pop a warning if the residue doesn't match, 
                    which usually means you are using a wrong residue id and requires a double check.)
                A : Chain index. This is optional witha default value: A.
                    (Different chains are defined by 'TER' marks in the PDB file and the index is noted in order.
                    * Note this will not always be the same as the chain id in the PDB line
                    * defined in self.stru)
                ##: Residue index. 
                    (Strictly correponding residue indexes in the original file.
                    * defined in self.stru)
                Y : Target residue letter.
                    (The residue type after the mutation.) 
                'WT'
                Can also use WT to specify the WT. The decoded MutaFlag will be ('WT', 'WT', 'WT', 'WT') to fit in the format.
            example:
                >>> pdb_obj.Add_MutaFlag('D83K')
                >>> pdb_obj.MutaFlags
                [('D', 'A', '83', 'K')]
                >>> pdb_obj.Add_MutaFlag('DA83K')
                >>> pdb_obj.MutaFlags
                [('D', 'A', '83', 'K')]
                >>> pdb_obj.Add_MutaFlag(['DA83K', 'EB226P'])
                >>> pdb_obj.MutaFlags
                [('D', 'A', '83', 'K'), ('E', 'B', '226', 'P')]
            
            **Assign random mutation(s)** 
                ['r', {'position': 'keyword', 
                       'target_resi':'keyword'}]
                To assigne random mutaion, start the list with 'r' or 'random' followed by a map that
                defines the rule of the randomlization. In the map dictionary there are 2 keys to fill:
                position        : use a keyword or a pattern to define availiable positions to mutate.
                target_residue  : use a keyword or a pattern to define availiable target residues.
                (* Note that when using random assignment, the list can no longer contain any manual assign str for mutation.
                e.g.: pdb_obj.Add_MutaFlag(['V23T', 'r']) is not valid)
            example:
                TODO(shaoqz): finish this function.
        if_U: 
            if include mutations to U (selenocysteine) in random generation.
        if_self:
            if "mutation to the same amino acid" is allowed for random mutation.

        Returns:
        There are two parts of actually outputs:
        In self.MutaFlags: Assigned mutations in a list of tuples
            example: 
            [('D', 'A', '83', 'K'), ('E', 'B', '226', 'P')]
        In return value: A str representing a tag of current mutations to the PDB.
            example:
            _DA83K_EB226P
        """

        if type(Flag) == str:

            if Flag == 'r' or Flag == 'random':
                resi_1 = ''
                resi_2 = ''
                Muta_c_id = ''
                Muta_r_id = ''

                # Flag Generation：
                # Random over the self.stru. Strictly correponding residue indexes in the original file. (no sort applied)
                self.get_stru()
                # random over the structure "protein" part.
                chain = choice(self.stru.chains)
                resi = choice(chain.residues)
                Muta_c_id = chain.id
                Muta_r_id = str(resi.id)
                if resi.name in Resi_map2:
                    resi_1 = Resi_map2[resi.name]
                else:
                    if Config.debug >= 1:
                        print('WARNING: pdb.Add_MutaFlag(): A non-canonical animo acid is being mutated! The MutaFlag will have a 3-letter code for the original residue .')
                    resi_1 = resi.name
                # random over the residue list
                if if_U:
                    m_Resi_list = Resi_list
                else:
                    m_Resi_list = Resi_list[:-1]
                resi_2 = choice(m_Resi_list)
                # avoid or not mutation to self
                if not if_self:
                    while resi_2 == resi_1:
                        resi_2 = choice(m_Resi_list)

                MutaFlag = (resi_1, Muta_c_id, Muta_r_id ,resi_2)
                self.MutaFlags.append(MutaFlag)

            else:
                MutaFlag = self._read_MutaFlag(Flag)
                self.MutaFlags.append(MutaFlag)

        if type(Flag) == list:
            for i in Flag:
                MutaFlag = self._read_MutaFlag(i)
                self.MutaFlags.append(MutaFlag)

        if Config.debug >= 1:
            print('Current MutaFlags:')
            for flag in self.MutaFlags:
                print(self._build_MutaName(flag))

        label=''
        for flag in self.MutaFlags:
            label=label+'_'+self._build_MutaName(flag) 
        return label
                

    def _read_MutaFlag(self, Flag):
        '''
        decode the manually input MutaFlag. Return (resi_1, chain_id, resi_id, resi_2)
        --------------
        Grammer(XA11Y):
        X : Original residue name. Leave X if unknow. 
            Only used for build filenames. **Do not affect any calculation.**
        A : Chain index. Determine by 'TER' marks in the PDB file. (Do not consider chain_indexs in the original file.)
        11: Residue index. Strictly correponding residue indexes in the original file. (NO sort applied)
        Y : Target residue name.

        * can also use WT to specify the WT. The decoded MutaFlag will be ('WT', 'WT', 'WT', 'WT') to fit in the format

        the later 3 will go through a san check to make sure they are within the range:
        A : chain.id range in self.stru.chains
        11: resi.id range in self.stru.chain[int].residues
        Y : Resi_list
        '''
        pattern = r'([A-Z])([A-Z])?([0-9]+)([A-Z])'
        if Flag == 'WT':
            return ('WT', 'WT', 'WT', 'WT')
        F_match = re.match(pattern, Flag)
        if F_match is None:
            raise Exception('_read_MutaFlag: Required format: XA123Y (or X123Y indicating the first chain)')
        
        resi_1 = F_match.group(1)
        chain_id = F_match.group(2)
        resi_id = str(F_match.group(3))
        resi_2 = F_match.group(4)

        # default
        if F_match.group(2) is None:
            chain_id = 'A'
            if Config.debug >= 1:
                print('_read_MutaFlag: No chain_id is provided! Mutate in the first chain by default. Input: ' + Flag)   

        # san check of the manual input
        self.get_stru()
        chain_id_list = [i.id for i in self.stru.chains]
        if not chain_id in chain_id_list:
            raise Exception('_read_MutaFlag: San check failed. Input chain id in not in range.'+line_feed+' range: '+ repr(chain_id_list))
        chain_int = ord(chain_id)-65
        resi_id_list = [str(i.id) for i in self.stru.chains[chain_int].residues]
        if not resi_id in resi_id_list:
            raise Exception('_read_MutaFlag: San check failed. Input resi id in not in range.'+line_feed+' range: '+ repr(resi_id_list))
        if not resi_2 in Resi_list:
            raise Exception('_read_MutaFlag: Only support mutate to the known 21 residues. AmberMaps.Resi_list: '+ repr(Resi_list))


        return (resi_1, chain_id, resi_id, resi_2)
    

    def _build_MutaName(self, Flag):
        '''
        Take a MutaFlag Tuple and return a str of name
        '''
        if 'WT' in Flag:
            return 'WT'
        return Flag[0]+Flag[1]+Flag[2]+Flag[3]


    def PDBMin(
        self,
        cycle: int = 20000,
        engine: str = 'Amber_GPU',
        if_cluster_job: bool = 1, 
        cluster: ClusterInterface = None, 
        period: int = 180,
        res_setting: Union[dict, None] = None,
        cpu_cores: Union[str, None] = None,
        cluster_debug: bool = 0
        ):
        '''
        Run a minization use self.prmtop and self.inpcrd and setting form class Config.
        --------------------------------------------------
        Save changed PDB to self.path (containing water and ions)
        Mainly used for remove bad contact from PDB2PDBwLeap.
        cpu_cores
            provide these for configure the gaussian input file when using a non-dict res_setting.
            these values should be the same as indicated in the res_setting.
        '''
        
        out4_PDB_path=self.path_name+'_min.pdb'

        #make MD config
        min_dir = self.cache_path+'/PDBMin'
        minin_path = min_dir + '/min.in'
        minout_path = min_dir + '/min.out'
        minrst_path = min_dir + '/min.ncrst' # change to ncrst seeking for solution of rst error
        mkdir(min_dir)
        self._build_MD_min(min_dir, cycle=cycle)

        # determine cluster rescources
        if if_cluster_job:
            core_type = engine.split('_')[-1] # this go front because needed in expressing cmd
            res_setting, env_settings = type(self)._prepare_cluster_job_pdbmd(cluster, core_type, res_setting)
            if core_type == 'CPU':
                # should some how wrap this up to a function
                if isinstance(res_setting, dict):
                    cpu_cores = res_setting['node_cores']
                else: # directly use argument -> support non-dict res_setting input e.g.: str
                    if cpu_cores == None:
                        raise TypeError('ERROR: Requires cpu_cores input for configuring sander if using CPU and provide res_setting not as a dict')
        elif cpu_cores != None:
            raise TypeError('ERROR: cpu_cores should be None if not submit to cluster.')
        # express engine
        PC_cmd, engine_path = Config.Amber.get_Amber_engine(engine=engine, n_cores=cpu_cores)

        # make cmd
        cmd = f'{PC_cmd} {engine_path} -O -i {minin_path} -o {minout_path} -p {self.prmtop_path} -c {self.inpcrd_path} -r {minrst_path} -ref {self.inpcrd_path}'
        
        # run
        if if_cluster_job:
            job = job_manager.ClusterJob.config_job(
                commands = cmd,
                cluster = cluster,
                env_settings = env_settings,
                res_keywords = res_setting,
                sub_dir = './', # because path are relative
                sub_script_path = f'{min_dir}/submit_PDBMin_{core_type}.cmd'
            )
            job.submit()
            if Config.debug > 0:
                print(f'''Running Min on {cluster.NAME}: job_id: {job.job_id} script: {job.sub_script_path} period: {period}''')
            job.wait_to_end(period=period)
            try:
                PDB._detect_amber_error(job)
            except AmberError as e:
                if not PDB._is_normal_amber_output(minout_path):
                    raise e
            
            if cluster_debug:
                self.pdbmin_job = job
        else:
            if Config.debug >= 1:
                print(f'running: {cmd}')
            os.system(cmd)
        
        # rst2pdb
        run('ambpdb -p '+self.prmtop_path+' -c '+minrst_path+' > '+out4_PDB_path, check=True, text=True, shell=True, capture_output=True)
        
        # clean
        os.system('mv '+self.prmtop_path+' '+self.inpcrd_path+' '+min_dir)

        # update
        self.path = out4_PDB_path
        self._update_name()
        self.prmtop_path=min_dir+'/'+self.prmtop_path.split('/')[-1]
        self.inpcrd_path=min_dir+'/'+self.inpcrd_path.split('/')[-1]
        return self.path


    @classmethod
    def _is_normal_amber_output(cls, out_file):
        '''check if the amber output file terminated normally'''
        with open(out_file) as f:
            if "5.  TIMINGS" in f.read():
                return True
            return False


    def rm_allH(self, ff='Amber', if_ligand=0):
        '''
        remove wrong hydrogens added by leap after mutation. (In the case that the input file was a H-less one from crystal.)
        ----------
        if_ligand:
        0 - remove Hs of standard protein residues only.
        1 - remove all Hs base on the nomenclature. (start with H and not in the non_H_list)
        '''
        # out path
        o_path=self.path_name+'_rmH.pdb'

        # crude judgement of H including customized H
        if if_ligand:
            not_H_list = ['HG', 'HF', 'HS'] # non-H elements that start with "H"
            with open(self.path) as f:
                with open(o_path,'w') as of:
                    for line in f:
                        if line.startswith('ATOM'):
                            atom_name = line[12:16].strip()
                            resi_name = line[17:20].strip()
                            if atom_name[0] == 'H':
                                if len(atom_name) >= 2:
                                    if atom_name[:2] not in not_H_list:
                                        continue
                                    elif resi_name != atom_name[:2]: # non-metal
                                        continue

                                else:
                                    continue
                        of.write(line)
        else:
            # H list (residue only)
            H_namelist=[]
            for name in Resi_Ele_map[ff]:
                if Resi_Ele_map[ff][name] == 'H':
                    H_namelist.append(name)

            with open(self.path) as f:
                with open(o_path,'w') as of:
                    for line in f:
                        if line[12:16].strip() in H_namelist and line[17:20].strip() in Resi_map2.keys():
                            continue
                        of.write(line)
        self.path=o_path
        self._update_name()       


    '''
    ========
    General MD
    ========
    '''

    def PDB2FF(self, 
        prm_out_path='', 
        o_dir='', 
        lig_method='AM1BCC', 
        renew_lig=0, 
        local_lig=1, 
        ifsavepdb=0, 
        igb=None, 
        ifsolve=1,
        if_prm_only=0,
        lig_net_charge_mapper=None
        ):
        '''
        PDB2FF(self, o_dir='')
        --------------------
        prm_out_path: output path of the prmtop file
        o_dir contral where the leap.in and leap.log go: has to contain a / at the end (e.g.: ./dir/)
        renew_lig: 0 use old ligand parm files if detected.
                   1 generate new ones. 
        local_lig: 0 export lig files to the workdir level in HTP jobs.
                   1 keep local 
        --------------------
        chains:
        ligand: - less junk files if your workflow contains a protonation step in advance.  
        metal:
        '''
        # check and generate self.stru
        self.get_stru()

        # build things seperately
        if local_lig:
            self.lig_dir = self.dir+'/ligands/'
            self.met_dir = self.dir+'/metalcenters/'
        else:
            self.lig_dir = self.dir+'/../ligands/'
            self.met_dir = self.dir+'/../metalcenters/'
        mkdir(self.lig_dir)
        mkdir(self.met_dir)

        # metalcenters_path = self.stru.build_metalcenters(self.met_dir)
        # parm
        ligand_parm_paths = self._ligand_parm(self.lig_dir, method=lig_method, renew=renew_lig, net_charge_mapper=lig_net_charge_mapper)
        # self._metal_parm(metalcenters_path)
        # combine
        if o_dir != '':
            mkdir(o_dir)
        self._combine_parm(ligand_parm_paths, prm_out_path=prm_out_path, o_dir=o_dir, ifsavepdb=ifsavepdb, igb=igb, ifsolve=ifsolve, if_prm_only=if_prm_only)
        if ifsavepdb:
            self.path = self.path_name+'_ff.pdb'
            self._update_name()
        
        return (self.prmtop_path,self.inpcrd_path)

    
    def _ligand_parm(self, lig_dir, method='AM1BCC', lig_charge_method='PYBEL', lig_charge_ph=7.0, renew=0, net_charge_mapper=None):
        '''
        Turn ligands to prepi (w/net charge), parameterize with parmchk
        return [(perpi_1, frcmod_1), ...]
        -----------
        method  : method use for ligand charge. Only support AM1BCC now.
        renew   : 0:(default) use old parm files if exist. 1: renew parm files everytime
        TODO check if the ligand is having correct name. (isolate the renaming function and also use in the class structure)
        * WARN: The parm file for ligand will always be like xxx/ligand_1.frcmod. Remember to enable renew when different object is sharing a same path.
        * BUG: Antechamber has a bug that if current dir has temp files from previous antechamber run (ANTECHAMBER_AC.AC, etc.) sqm will fail. Now remove them everytime.
        '''
        parm_paths = []
        self.prepi_path = {}
        self.frcmod_path = {}
        
        lig_list = self.stru.get_all_ligands(ifunique=1)
        for lig in lig_list:
            lig: Ligand
            if Config.debug >= 1:
                print( 'Working on: '+'_'.join((lig.name, str(lig.id))) )            
            # target files
            out_prepi = lig_dir+'ligand_'+lig.name+'.prepin'
            out_frcmod = lig_dir+'ligand_'+lig.name+'.frcmod'
            # if renew
            if os.path.isfile(out_prepi) and os.path.isfile(out_frcmod) and not renew:
                if Config.debug >= 1:
                    print('Parm files exist: ' + out_prepi + ' ' + out_frcmod)
                    print('Using old parm files.')
            else:
                # build ligand pdb file
                lig_pdb_path = lig_dir+'ligand_'+lig.name+'.pdb'
                lig.build(lig_pdb_path, ft='PDB')
                # get net charge
                if not (net_charge_mapper and net_charge_mapper.get(lig.name, None) != None):
                    net_charge = lig.get_net_charge(method=lig_charge_method, ph=lig_charge_ph, o_dir=lig_dir)
                else:
                    net_charge = net_charge_mapper[lig.name]
                # get parameters
                if method == 'AM1BCC':
                    #gen prepi (net charge and correct protonation state is important)
                    if Config.debug >= 1:
                        print('running: '+Config.Amber.AmberHome+'/bin/antechamber -i '+lig_pdb_path+' -fi pdb -o '+out_prepi+' -fo prepi -c bcc -s 0 -nc '+str(net_charge))
                    run(Config.Amber.AmberHome+'/bin/antechamber -i '+lig_pdb_path+' -fi pdb -o '+out_prepi+' -fo prepi -c bcc -s 0 -nc '+str(net_charge), check=True, text=True, shell=True, capture_output=True)
                    if Config.debug <= 1:
                        os.system('rm ANTECHAMBER* ATOMTYPE.INF NEWPDB.PDB PREP.INF sqm.pdb sqm.in sqm.out')
                    #gen frcmod
                    if Config.debug >= 1:
                        print('running: '+Config.Amber.AmberHome+'/bin/parmchk2 -i '+out_prepi+' -f prepi -o '+out_frcmod)
                    run(Config.Amber.AmberHome+'/bin/parmchk2 -i '+out_prepi+' -f prepi -o '+out_frcmod, check=True, text=True, shell=True, capture_output=True)                
            #record
            parm_paths.append((out_prepi, out_frcmod))
            self.prepi_path[lig.name] = out_prepi
            self.frcmod_path[lig.name] = out_frcmod

        return parm_paths


    def _combine_parm(self, lig_parms, prm_out_path='', o_dir='', ifsavepdb=0, ifsolve=1, box_type=None, box_size=Config.Amber.box_size, igb=None, if_prm_only=0):
        '''
        combine different parmeter files and make finally inpcrd and prmtop
        -------
        structure: pdb
        ligands: prepi, frcmod
        metalcenters, artificial residues: TODO
        '''
        if box_type == None:
            box_type = Config.Amber.box_type
            
        leap_path= self.cache_path+'/leap.in'
        sol_path= self.path_name+'_ff.pdb'
        with open(leap_path, 'w') as of:
            of.write('source leaprc.protein.ff14SB'+line_feed)
            of.write('source leaprc.gaff'+line_feed)
            of.write('source leaprc.water.tip3p'+line_feed)
            # ligands
            for prepi, frcmod in lig_parms:
                of.write('loadAmberParams '+frcmod+line_feed)
                of.write('loadAmberPrep '+prepi+line_feed)
            of.write('a = loadpdb '+self.path+line_feed)
            if self.disulfied_residue_pairs:
                for ss_bond_pairs in self.disulfied_residue_pairs:
                    of.write(f'bond a.{ss_bond_pairs[0][1]}.SG a.{ss_bond_pairs[1][1]}.SG{line_feed}')
            # igb Radii
            if igb != None:
                radii = radii_map[str(igb)]
                of.write('set default PBRadii '+ radii +line_feed)
            of.write('center a'+line_feed)
            # solvation
            if ifsolve:
                of.write('addions a Na+ 0'+line_feed)
                of.write('addions a Cl- 0'+line_feed)
                if box_type == 'box':
                    of.write('solvatebox a TIP3PBOX '+box_size+line_feed)
                if box_type == 'oct':
                    of.write('solvateOct a TIP3PBOX '+box_size+line_feed)
                if box_type != 'box' and box_type != 'oct':
                    raise Exception('PDB._combine_parm().box_type: Only support box and oct now!')
            # save
            if prm_out_path == '':
                if o_dir == '':                        
                    of.write('saveamberparm a '+self.path_name+'.prmtop '+self.path_name+'.inpcrd'+line_feed)
                    self.prmtop_path=self.path_name+'.prmtop'
                    self.inpcrd_path=self.path_name+'.inpcrd'
                else:
                    of.write('saveamberparm a '+o_dir+self.name+'.prmtop '+o_dir+self.name+'.inpcrd'+line_feed)
                    self.prmtop_path=o_dir+self.name+'.prmtop'
                    self.inpcrd_path=o_dir+self.name+'.inpcrd'
            else:
                if o_dir == '':
                    if if_prm_only:
                        temp_dir = f'{self.dir}/temp'
                        mkdir(temp_dir)
                        of.write('saveamberparm a '+prm_out_path+f' {temp_dir}/tmp.inpcrd'+line_feed)
                        self.prmtop_path=prm_out_path
                        self.inpcrd_path=None
                    else:
                        of.write('saveamberparm a '+prm_out_path+' '+self.path_name+'.inpcrd'+line_feed)
                        self.prmtop_path=prm_out_path
                        self.inpcrd_path=self.path_name+'.inpcrd'
                else:
                    if if_prm_only:
                        temp_dir = f'{self.dir}/temp'
                        mkdir(temp_dir)
                        of.write('saveamberparm a '+prm_out_path+f' {temp_dir}/tmp.inpcrd'+line_feed)
                        self.prmtop_path=prm_out_path
                        self.inpcrd_path=None
                    else:
                        of.write('saveamberparm a '+prm_out_path+' '+o_dir+self.name+'.inpcrd'+line_feed)
                        self.prmtop_path=prm_out_path
                        self.inpcrd_path=o_dir+self.name+'.inpcrd'

            if ifsavepdb:
                of.write('savepdb a '+sol_path+line_feed)
            of.write('quit'+line_feed)

        try:
            run('tleap -s -f '+leap_path+' > '+leap_path[:-2]+'out', check=True,  text=True, shell=True, capture_output=True)
        except SubprocessError as e:
            print(f'stderr: {str(e.stderr).strip()}')
            print(f'stdout: {str(e.stdout).strip()}')
            raise e

        return self.prmtop_path, self.inpcrd_path


    def rm_wat(self):
        '''
        Remove water and ion for the pdb. Remians the same if there's no water or ion.
        Now only skip [Na+,Cl-,WAT,HOH] // Append more in the future.
        Save changed files into self.path.
        TODO: need to support key water.
        '''
        out_path = self.path_name+'_rmW.pdb'
        self._get_file_path()
        with open(self.path) as f:
            with open(out_path,'w') as of:
                skip_list=['Na+','Cl-','WAT','HOH']
                change_flag=0
                for line in f:
                    PDB_l = PDB_line(line)
                    #Skip some lines
                    skip_flag=0
                    #skip the CRYST1 line
                    if line[:6] == 'CRYST1':
                        change_flag=1
                        continue
                    #keep TER and END
                    if line[:3] == 'END':
                        of.write(line)
                        continue
                    if line[:3] == 'TER':
                        # if ter_flag is still 1. It means this first line kept after last TER is still a TER
                        # skip TER in this case
                        if ter_flag == 1:
                            continue
                        of.write(line)
                        ter_flag = 1
                        continue

                    #skip the water and ion(Append in the future)
                    for i in skip_list:
                        if PDB_l.resi_name == i:
                            skip_flag=1
                            break
                    if skip_flag:
                        change_flag=1
                        continue

                    of.write(line)
                    # set to 0 until first following line after TER is kept (not TER or skip_list)
                    ter_flag=0

                if not change_flag:
                    print('rm_wat(): No change.')

        self.path=out_path
        self._update_name()

        return self.path


    def PDBMD(  
        self, 
        tag: str = '', 
        o_dir: str = '', 
        engine: str = 'Amber_GPU', 
        equi_cpu: bool = 0,
        if_cluster_job: bool = 1, # TODO(qz) maybe its better to make everything below a dictionary 
        cluster: ClusterInterface = None, # since the requirement is similar for a cluster job
        period: int = 600,
        res_setting: Union[dict, None] = None,
        cpu_cores: Union[str, int, None] = None,
        res_setting_equi_cpu: Union[dict, None] = None,
        equi_cpu_cores: Union[str, int, None] = None,
        cluster_debug: bool = 0,
        check_out_put: bool = 1,
        ):
        '''
        Use self.prmtop_path and self.inpcrd_path to initilize a MD simulation.
        The default MD configuration settings are assigned by class Config.Amber.
        * User can also set MD configuration for the current object by assigning values in self.conf_xxxx.
        * e.g.: self.conf_heat['nstlim'] = 50000
        --------------
        Args:
        self    :
            dir 
                - for creating the MD storage dir under
            prmtop_path 
                - MD input file
            inpcrd_path 
                - MD input file
        o_dir   : 
            Write files in o_dir (current self.dir/MD by default).
        tag     : 
            tag the assigned tag to the MD folder
        engine  : 
            MD engine (cpu/gpu)
        equi_cpu: 
            if use cpu for equi step
        if_cluster_job:
            if submit the calculation to a cluster using the job manager
        cluster: (if_cluster_job is 1)
            The cluster used
        period: (if_cluster_job is 1)
            the time cycle for each job state check (Unit: s) (default: 600)
        res_setting: (if_cluster_job is 1)
            (default: Config.Amber.MD_CPU_RES)
            provide specific keys you want to change the rest will remain default.
            (e.g. res_setting = {'node_cores':'2'} will still use a complete dictionary 
            with only node_cores modified.)
        res_setting_equi_cpu: 
            resource for equi cpu job if use cpu for equi
            provide specific keys you want to change the rest will remain default.
        cpu_cores, equi_cpu_cores
            provide these for configure the gaussian input file when using a non-dict res_setting.
            these values should be the same as indicated in the res_setting.
    cluster_debug:
           1:  add also the pdbmd job obj to the pdb obj
        --data--
        Attribute:
            store resulting NetCrd traj file path in self.nc
        Return:
            Return the nc path of the prod step
        '''
        # make folder
        if o_dir == '':
            o_dir = self.dir+'/MD'+tag
        mkdir(o_dir)
        # default cpu cores
        # format settings for cluster job (since it also influence the command)
        if if_cluster_job:
            core_type = engine.split('_')[-1]
            res_setting, env_settings = type(self)._prepare_cluster_job_pdbmd(cluster, core_type, res_setting)
            if core_type == 'CPU':
                # TODO(shaoqz) see PDBMin, somehow wrap this up
                if isinstance(res_setting, dict):
                    cpu_cores = res_setting['node_cores']
                else: # support non-dict res_setting input e.g.: str
                    if cpu_cores == None:
                        raise TypeError('ERROR: Requires cpu_cores input for configuring sander if using CPU and provide res_setting not as a dict')

            if core_type == 'GPU' and equi_cpu:
                # get env res for equi
                env_settings_equi_cpu = cluster.AMBER_ENV['CPU']
                res_setting_equi_cpu = type(self)._get_default_res_setting_pdbmd(res_setting_equi_cpu, 'CPU')
                if isinstance(res_setting_equi_cpu, dict):
                    equi_cpu_cores = res_setting_equi_cpu['node_cores']
                else: # support non-dict res_setting input e.g.: str
                    if equi_cpu_cores == None:
                        raise TypeError('ERROR: Requires equi_cpu_cores input for configuring sander if using CPU and provide equi_res_setting not as a dict')
        else:
            if cpu_cores != None or equi_cpu_cores != None:
                raise TypeError('ERROR: cpu_cores or equi_cpu_cores should be None if not submit to cluster.')        

        # express engine (pirority: AmberEXE_GPU/AmberEXE_CPU - AmberHome/bin/xxx)
        PC_cmd, engine_path = Config.Amber.get_Amber_engine(engine=engine, n_cores=cpu_cores)
        # express cpu engine if equi_cpu
        if equi_cpu:
            if Config.Amber.AmberEXE_CPU == None:
                cpu_engine_path = Config.Amber.AmberHome+'/bin/sander.MPI'
            else:
                cpu_engine_path = Config.Amber.AmberEXE_CPU

        # build input file (use self.MD_conf_xxxx)
        min_path = self._build_MD_min(o_dir)
        heat_path = self._build_MD_heat(o_dir)
        equi_path = self._build_MD_equi(o_dir)
        prod_path = self._build_MD_prod(o_dir)
        # build md commands
        cmd_min = f'{PC_cmd} {engine_path} -O -i {min_path} -o {o_dir}/min.out -p {self.prmtop_path} -c {self.inpcrd_path} -r {o_dir}/min.rst -ref {self.inpcrd_path}{line_feed}'.strip(' ')
        cmd_heat = f'{PC_cmd} {engine_path} -O -i {heat_path} -o {o_dir}/heat.out -p {self.prmtop_path} -c {o_dir}/min.rst -r {o_dir}/heat.rst -ref {o_dir}/min.rst{line_feed}'.strip(' ')
        cmd_prod = f'{PC_cmd} {engine_path} -O -i {prod_path} -o {o_dir}/prod.out -p {self.prmtop_path} -c {o_dir}/equi.rst -r {o_dir}/prod.rst -ref {o_dir}/equi.rst -x {o_dir}/prod.nc{line_feed}'.strip(' ')
        # gpu debug for equi
        if equi_cpu:
            cmd_equi = f'{Config.get_PC_cmd(equi_cpu_cores)} {cpu_engine_path} -O -i {equi_path} -o {o_dir}/equi.out -p {self.prmtop_path} -c {o_dir}/heat.rst -r {o_dir}/equi.rst -ref {o_dir}/heat.rst -x {o_dir}/equi.nc{line_feed}'.strip(' ')
        else:
            cmd_equi = f'{PC_cmd} {engine_path} -O -i {equi_path} -o {o_dir}/equi.out -p {self.prmtop_path} -c {o_dir}/heat.rst -r {o_dir}/equi.rst -ref {o_dir}/heat.rst -x {o_dir}/equi.nc{line_feed}'.strip(' ')

        # run MD exe
        if if_cluster_job:
            md_jobs = []
            if equi_cpu:
                # divide into 3 jobs
                cmd_1 = cmd_min + cmd_heat
                cmd_2 = cmd_equi
                cmd_3 = cmd_prod
                # make job: 1&3 are in same environment 2 is always in CPU environment
                job_1 = job_manager.ClusterJob.config_job(
                    commands = cmd_1,
                    cluster = cluster,
                    env_settings = env_settings,
                    res_keywords = res_setting,
                    sub_dir = './', # because path are relative
                    sub_script_path = f'{o_dir}/submit_PDBMD_1_{core_type}.cmd')
                job_1.submit()
                if Config.debug > 0:
                    print(f'''Running MD on {cluster.NAME}: job_id: {job_1.job_id} script: {job_1.sub_script_path} period: {period}''')
                job_1.wait_to_end(period)
                type(self)._detect_amber_error(job_1)

                job_2 = job_manager.ClusterJob.config_job(
                    commands = cmd_2,
                    cluster = cluster,
                    env_settings = env_settings_equi_cpu,
                    res_keywords = res_setting_equi_cpu,
                    sub_dir = './', # because path are relative
                    sub_script_path = f'{o_dir}/submit_PDBMD_2_CPU.cmd')
                job_2.submit()
                if Config.debug > 0:
                    print(f'''Running MD on {cluster.NAME}: job_id: {job_2.job_id} script: {job_2.sub_script_path} period: {period}''')
                job_2.wait_to_end(period)
                type(self)._detect_amber_error(job_2)

                job_3 = job_manager.ClusterJob.config_job(
                    commands = cmd_3,
                    cluster = cluster,
                    env_settings = env_settings,
                    res_keywords = res_setting,
                    sub_dir = './', # because path are relative
                    sub_script_path = f'{o_dir}/submit_PDBMD_3_{core_type}.cmd')
                job_3.submit()
                if Config.debug > 0:
                    print(f'''Running MD on {cluster.NAME}: job_id: {job_3.job_id} script: {job_3.sub_script_path} period: {period}''')
                job_3.wait_to_end(period)
                type(self)._detect_amber_error(job_3)
                md_jobs.extend([job_1, job_2, job_3])
            else:
                # all GPU or CPU
                cmd = cmd_min + cmd_heat + cmd_equi + cmd_prod
                job = job_manager.ClusterJob.config_job(
                    commands = cmd,
                    cluster = cluster,
                    env_settings = env_settings,
                    res_keywords = res_setting,
                    sub_dir = './', # because path are relative
                    sub_script_path = f'{o_dir}/submit_PDBMD_{core_type}.cmd'
                )
                job.submit()
                if Config.debug > 0:
                    print(f'''Running MD on {cluster.NAME}: job_id: {job.job_id} script: {job.sub_script_path} period: {period}''')
                job.wait_to_end(period=period)
                type(self)._detect_amber_error(job)
                md_jobs.append(job)
            
            if cluster_debug:
                self.pdbmd_jobs = md_jobs

        else:
            if Config.debug >= 1:
                print(f'running: {cmd_min}')
            os.system(cmd_min)
            if Config.debug >= 1:
                print(f'running: {cmd_heat}')
            os.system(cmd_heat)
            if Config.debug >= 1:
                print(f'running: {cmd_equi}')
            os.system(cmd_equi)
            if Config.debug >= 1:
                print(f'running: {cmd_prod}')
            os.system(cmd_prod)
        if check_out_put and not PDB._is_normal_amber_output(f'{o_dir}/prod.out'):
            raise Exception("Amber production did not terminate normally")

        # return value
        self.nc = o_dir+'/prod.nc'
        return o_dir+'/prod.nc'

    @staticmethod
    def _detect_amber_error(amber_job):
        '''
        amber pmemd tend to delay the error show up in the workflow.
        It does not provide abnormal exit code but just output some error information to the stdout/stderr
        '''
        with open(amber_job.job_cluster_log) as f:
            f_str = f.read()
            if 'ERROR' in f_str or 'Error' in f_str:
                raise AmberError(f'Amber Terminated Abnormally! Here is the output ({amber_job.job_cluster_log}):{line_feed} {f_str}') 
                # TODO make a workflow exception that can show some more step releted info
                # So that can use try except to do some special action when dont want one fail ruin the HTP
                # this may also catch unnessessay error like MPI ones

    @staticmethod
    def _get_default_res_setting_pdbmd(res_setting, core_type):
        if res_setting is None:
            res_setting = dict()
        if isinstance(res_setting, dict):
            res_setting_holder = copy.deepcopy(Config.Amber.MD_RES[core_type])
            # replace assigned key in default dict
            for k, v in res_setting.items():
                res_setting_holder[k] = v
            return res_setting_holder
        if isinstance(res_setting, str):
            return res_setting

    @classmethod
    def _prepare_cluster_job_pdbmd(cls, cluster, core_type, res_setting) -> tuple[str, str]:
        '''
        update the default res_setting with the input
        get env_setting according to the core type
        '''
        #san check
        if not isinstance(cluster, ClusterInterface):
            raise TypeError('cluster job need a cluster (ClusterInterface object) input')
        # get res and env settinng
        # interface check
        if 'AMBER_ENV' not in dir(cluster):
            raise Exception('PDBMD requires the input cluster have the AMBER_ENV attr')
        env_settings = cluster.AMBER_ENV[core_type]
        # default for res
        res_setting = cls._get_default_res_setting_pdbmd(res_setting, core_type)

        return res_setting, env_settings


    def _build_MD_min(self, o_dir, cycle=None):
        '''
        Build configuration file for a minimization job
        See default value Config.Amber.conf_min
        '''
        #path
        o_path=o_dir+'/min.in'
        #maxcyc related
        if cycle:
            maxcyc = cycle
        else:
            maxcyc = self.conf_min['maxcyc']
        if self.conf_min['ncyc'] == '0.5maxcyc':
            ncyc = str(int(0.5 * maxcyc))
        if self.conf_min['ntpr'] == '0.01maxcyc':
            ntpr = str(int(0.01 * maxcyc))
        maxcyc = str(maxcyc)
        #restrain related
        if self.conf_min['ntr'] == '1':
            ntr_line = '  ntr   = '+self.conf_min['ntr']+',	 restraint_wt = '+self.conf_min['restraint_wt']+', restraintmask = '+self.conf_min['restraintmask']+','+line_feed
        else:
            ntr_line = ''
        if self.conf_min['nmropt_rest'] == '1':
            self.conf_min['nmropt'] = '1'
            nmropt_line = '  nmropt= '+self.conf_min['nmropt']+','+line_feed
            DISANG_tail = ''' &wt
  type='END'
 /
  DISANG= '''+self.conf_min['DISANG']+line_feed
            self._build_MD_rs(step='min',o_path=self.conf_min['DISANG'])
        else:
            nmropt_line = ''
            DISANG_tail = ''
        #text        
        conf_str='''Minimize
 &cntrl
  imin  = 1,  ntx   = 1,  irest = 0,
  ntc   = '''+self.conf_min['ntc']+''',    ntf = '''+self.conf_min['ntf']+''',
  cut   = '''+self.conf_min['cut']+''',
  maxcyc= '''+maxcyc+''', ncyc  = '''+ncyc+''',
  ntpr  = '''+ntpr+''', ntwx  = 0,
'''+ntr_line+nmropt_line+''' /
'''+DISANG_tail
        #write
        with open(o_path,'w') as of:
            of.write(conf_str)
        return o_path


    def _build_MD_heat(self, o_dir):
        '''
        Build configuration file for a heat job
        See default value Config.Amber.conf_heat
        '''
        #path
        o_path=o_dir+'/heat.in'

        # nstlim related
        nstlim=self.conf_heat['nstlim']
        if self.conf_heat['A_istep2'] == '0.9nstlim':
            A_istep2=str(int(nstlim*0.9))
        if self.conf_heat['B_istep1'] == 'A_istep2+1':
            B_istep1=str(int(A_istep2)+1)
        if self.conf_heat['ntpr'] == '0.01nstlim':
            ntpr = str(int(nstlim*0.01))
        if self.conf_heat['ntwx'] == 'nstlim':
            ntwx = str(nstlim)
        nstlim = str(nstlim)
        #restrain related
        if self.conf_heat['ntr'] == '1':
            ntr_line = '  ntr   = '+self.conf_heat['ntr']+', restraint_wt = '+self.conf_heat['restraint_wt']+', restraintmask = '+self.conf_heat['restraintmask']+','+line_feed
        else:
            ntr_line = ''
        if self.conf_heat['nmropt_rest'] == '1':
            DISANG_tail = '''  DISANG='''+self.conf_heat['DISANG']+line_feed
            self._build_MD_rs(step='heat',o_path=self.conf_heat['DISANG'])
        else:
            DISANG_tail = ''
        conf_str='''Heat
 &cntrl
  imin  = 0,  ntx = 1, irest = 0,
  ntc   = '''+self.conf_heat['ntc']+''', ntf = '''+self.conf_heat['ntf']+''',
  cut   = '''+self.conf_heat['cut']+''',
  nstlim= '''+nstlim+''', dt= '''+self.conf_heat['dt']+''',
  tempi = '''+self.conf_heat['tempi']+''',  temp0='''+self.conf_heat['temp0']+''',  
  ntpr  = '''+ntpr+''',  ntwx='''+ntwx+''',
  ntt   = '''+self.conf_heat['ntt']+''', gamma_ln = '''+self.conf_heat['gamma_ln']+''',
  ntb   = '''+self.conf_heat['ntb']+''',  ntp = '''+self.conf_heat['ntp']+''',
  iwrap = '''+self.conf_heat['iwarp']+''',
  nmropt= 1,
  ig    = -1,
'''+ntr_line+''' /
 &wt
  type  = 'TEMP0',
  istep1= 0, istep2='''+A_istep2+''',
  value1= '''+self.conf_heat['tempi']+''', value2='''+self.conf_heat['temp0']+''',
 /
 &wt
  type  = 'TEMP0',
  istep1= '''+B_istep1+''', istep2='''+nstlim+''',
  value1= '''+self.conf_heat['temp0']+''', value2='''+self.conf_heat['temp0']+''',
 /
 &wt
  type  = 'END',
 /
'''+DISANG_tail
        #write
        with open(o_path,'w') as of:
            of.write(conf_str)
        return o_path


    def _build_MD_equi(self, o_dir):
        '''
        Build configuration file for the equilibration step
        See default value Config.Amber.conf_equi
        default ntwx -> 10ps
        '''
        #path
        o_path=o_dir+'/equi.in'

        # nstlim related
        nstlim=self.conf_equi['nstlim']
        if self.conf_equi['ntpr'] == '0.002nstlim':
            ntpr = str(int(nstlim*0.002))
        nstlim = str(nstlim)
        #restrain related
        if self.conf_equi['ntr'] == '1':
            ntr_line = '  ntr   = '+self.conf_equi['ntr']+', restraint_wt = '+self.conf_equi['restraint_wt']+', restraintmask = '+self.conf_equi['restraintmask']+','+line_feed
        else:
            ntr_line = ''
        if self.conf_equi['nmropt_rest'] == '1':
            self.conf_equi['nmropt'] = '1'
            nmropt_line = '  nmropt= '+self.conf_equi['nmropt']+','+line_feed
            DISANG_tail = ''' &wt
  type='END'
 /
  DISANG= '''+self.conf_equi['DISANG']+line_feed
            self._build_MD_rs(step='equi',o_path=self.conf_equi['DISANG'])
        else:
            nmropt_line = ''
            DISANG_tail = ''


        conf_str='''Equilibration:constant pressure
 &cntrl
  imin  = 0,  ntx = '''+self.conf_equi['ntx']+''',  irest = '''+self.conf_equi['irest']+''',
  ntf   = '''+self.conf_equi['ntf']+''',  ntc = '''+self.conf_equi['ntc']+''',
  nstlim= '''+nstlim+''', dt= '''+self.conf_equi['dt']+''',
  cut   = '''+self.conf_equi['cut']+''',
  temp0 = '''+self.conf_equi['temp0']+''',
  ntpr  = '''+ntpr+''', ntwx = '''+self.conf_equi['ntwx']+''',
  ntt   = '''+self.conf_equi['ntt']+''', gamma_ln = '''+self.conf_equi['gamma_ln']+''',
  ntb   = '''+self.conf_equi['ntb']+''',  ntp = '''+self.conf_equi['ntp']+''',
  iwrap = '''+self.conf_equi['iwarp']+''',
  ig    = -1,
'''+ntr_line+nmropt_line+''' /
'''+DISANG_tail
        #write
        with open(o_path,'w') as of:
            of.write(conf_str)
        return o_path


    def _build_MD_prod(self, o_dir):
        '''
        Build configuration file for the production step
        See default value Config.Amber.conf_prod
        default ntwx -> 10ps
        '''
        #path
        o_path=o_dir+'/prod.in'

        # nstlim related
        nstlim=self.conf_prod['nstlim']
        if self.conf_prod['ntpr'] == '0.001nstlim':
            ntpr = str(int(nstlim*0.001))
        nstlim = str(nstlim)
        #restrain related
        if self.conf_prod['ntr'] == '1':
            ntr_line = '  ntr   = '+self.conf_prod['ntr']+', restraint_wt = '+self.conf_prod['restraint_wt']+', restraintmask = '+self.conf_prod['restraintmask']+','+line_feed
        else:
            ntr_line = ''
        if self.conf_prod['nmropt_rest'] == '1':
            self.conf_prod['nmropt'] = '1'
            nmropt_line = '  nmropt= '+self.conf_prod['nmropt']+','+line_feed
            DISANG_tail = ''' &wt
  type='END'
 /
  DISANG= '''+self.conf_prod['DISANG']+line_feed
            self._build_MD_rs(step='prod',o_path=self.conf_prod['DISANG'])
        else:
            nmropt_line = ''
            DISANG_tail = ''

        conf_str='''Production
 &cntrl
  imin  = 0, ntx = '''+self.conf_prod['ntx']+''', irest = '''+self.conf_prod['irest']+''',
  ntf   = '''+self.conf_prod['ntf']+''',  ntc = '''+self.conf_prod['ntc']+''',
  nstlim= '''+nstlim+''', dt= '''+self.conf_prod['dt']+''',
  cut   = '''+self.conf_prod['cut']+''',
  temp0 = '''+self.conf_prod['temp0']+''',
  ntpr  = '''+ntpr+''', ntwx = '''+self.conf_prod['ntwx']+''',
  ntt   = '''+self.conf_prod['ntt']+''', gamma_ln = '''+self.conf_prod['gamma_ln']+''',
  ntb   = '''+self.conf_prod['ntb']+''',  ntp = '''+self.conf_prod['ntp']+''',
  iwrap = '''+self.conf_prod['iwarp']+''',
  ig    = -1,
'''+ntr_line+nmropt_line+''' /
'''+DISANG_tail
        #write
        with open(o_path,'w') as of:
            of.write(conf_str)
        return o_path
    

    def _build_MD_rs(self,step,o_path):
        '''
        Generate a file for DISANG restraint. Get parameters from self.conf_step.
        '''
        rs_str=''
        rs_data_step=self.__dict__['conf_'+step]['rs_constraints']
        for rest_data in rs_data_step:
            rs_str=rs_str+'''  &rst
   iat=  '''+','.join(rest_data['iat'])+''', r1= '''+rest_data['r1']+''', r2= '''+rest_data['r2']+''', r3= '''+rest_data['r3']+''', r4= '''+rest_data['r4']+''',
   rk2='''+rest_data['rk2']+''', rk3='''+rest_data['rk3']+''', ir6='''+rest_data['ir6']+''', ialtd='''+rest_data['ialtd']+''',
  &end
'''
        #write
        with open(o_path,'w') as of:
            of.write(rs_str)
        return o_path


    def reset_MD_conf(self):
        '''
        reset MD configuration of current object to default. 
        '''
        self.conf_min = Config.Amber.conf_min
        self.conf_heat = Config.Amber.conf_heat
        self.conf_equi = Config.Amber.conf_equi
        self.conf_prod = Config.Amber.conf_prod


    def show_MD_conf(self):
        '''
        Show MD configuration of current object. 
        '''
        print('Min :     '+repr(self.conf_min))
        print('Heat:     '+repr(self.conf_heat))
        print('Equi:     '+repr(self.conf_equi))
        print('Prod:     '+repr(self.conf_prod))

    '''
    ========
    QM/MM
    ========
    '''
    def PDB2QMMM(self, o_dir='',tag='', work_type='spe', qm='g16', keywords='', prmtop_path=None, prepi_path:dict=None, spin_list=[1,1], ifchk=1):
        '''
        generate QMMM input template based on [connectivity, atom order, work type, layer/freeze settings, charge settings]
        * NEED TO SET LAYER BY ATOM INDEX OR SELECT A LAYER PRESET (use Config.Gaussian.layer_preset and Config.Gaussian.layer_atoms)
        * define self.frames in the func
        --------
        qm          : QM program (default: g16 / Gaussian16)
        work_type   : QMMM calculation type (default: spe)
        o_dir       : out put directory of the input file (default: self.dir/QMMM/ )
        tag         : folder name tag for potential multiple mutations
        keywords    : additional keywords add to the work_type correlated route
        prmtop_path : provide prmtop file for determining charge and spin. use self.prmtop by default
        prepi_path  : a diction of prepin file path with each ligand name as key. (e.g.: {'4CO':'./ligand/xxx.prepin'})
        spin_list   : a list of spin for each layers. (Do not support auto judge of the spin now)
        ifchk       : if save chk and return chk paths
        <see more options in Config.Gaussian>
        ========
        Gaussian
        ========
        # route section (from config module / leave work type as an inp arg)
        # charge and spin
        # coordinate 
        	- atom label (from .lib)
        	- atom charge
        	- freeze part (general option / some presupposition / freeze MM) 
        	- layer (same as above)
        	- xyz (1. new system 2. existing template (do we still need?) -- from pdb / mdcrd / gout) ref: ONIOM_template_tool
        # connectivity
        # missing parameters (ligand related?)
        ---------
        In Config.Gaussian 
        ---------
        n_cores     : Cores for gaussian job (higher pirority)
        max_core    : Per core memory in MB for gaussian job (higher pirority)
        keywords    : a list of keywords that joined together when build
        layer_preset: default preset id copied to self.layer_preset
        layer_atoms : default layer_atoms copied to self.layer_atoms
        
        *om_lvl     : *(can only be edit manually before loading the module) oniom method level
        '''
        # san check
        support_work_type=['spe','opt','tsopt']
        if work_type not in support_work_type:
            raise Exception('PDB2QMMM.work_type : only support: '+repr(support_work_type))
        support_qm=['g16']
        if qm not in support_qm:
            raise Exception('PDB2QMMM.qm: only support: '+repr(support_qm))
        # default
        if prmtop_path == None:
            prmtop_path = self.prmtop_path
        # make folder
        if o_dir == '':
            o_dir = self.dir+'/QMMM'+tag
        mkdir(o_dir)
        # file path and name
        o_name = self.name+'_QMMM'
        g_temp_path = o_dir+'/'+o_name+'.gjf'
        #get stru
        self.get_stru()
        #get layer
        self._get_oniom_layer()
        # prepin path
        if prepi_path == None:
            prepi_path = self.prepi_path

        # build template
        if qm == 'g16':
            self.route = self._get_oniom_g16_route(work_type, o_name, key_words=keywords)
            title = 'ONIOM input template generated by PDB2QMMM module of XXX(software name)'+line_feed
            chrgspin = self._get_oniom_chrgspin(prmtop_path=prmtop_path, spin_list=spin_list)
            cnt_table = self.stru.get_connectivty_table(prepi_path=prepi_path)
            coord = self._get_oniom_g16_coord(prmtop_path) # use connectivity info from the line above.
            add_prm = self._get_oniom_g16_add_prm() # test for rules of missing parameters
            
            #combine and write 
            with open(g_temp_path,'w') as of:
                of.write(self.route) 
                of.write(line_feed) 
                of.write(title)
                of.write(line_feed)
                of.write(chrgspin)
                of.write(line_feed)
                of.write(coord)
                of.write(line_feed)
                of.write(cnt_table)
                of.write(line_feed)
                of.write(add_prm)
        
        # deploy to inp files
        frames = Frame.fromMDCrd(self.mdcrd)
        self.frames = frames
        gjf_paths = []
        chk_paths = []
        if Config.debug >= 1:
            print('Writing QMMM gjfs.')
        for i, frame in enumerate(frames):
            if ifchk:
                frame_path = frame.write_to_template(g_temp_path, index = str(i), ifchk=1)
                gjf_paths.append(frame_path[0])
                chk_paths.append(frame_path[1])
            else:
                gjf_paths.append(frame.write_to_template(g_temp_path, index = str(i), ifchk=0))
        # run Gaussian job
        self.qmmm_out = PDB.Run_QM(gjf_paths)

        if ifchk:
            self.qmmm_chk = chk_paths
            return self.qmmm_out, self.qmmm_chk

        return self.qmmm_out


    def _get_oniom_layer(self):
        '''
        get oniom layer base on self.layer_atoms or self.layer_preset
        save a Layer object to self.layer
        '''
        #  san check (need at least a set or a preset mode)
        if self.layer_atoms == [] and self.layer_preset == 0:
            raise Exception('PDB2QMMM: need layer setting. Please use self.set_oniom_layer.')
        # layer_atoms are in higher pirority
        if self.layer_atoms != []:
            self.layer = Layer(self, self.layer_atoms)
        else:
            self.layer = Layer.preset(self, self.layer_preset)


    def _get_oniom_g16_route(self, work_type, chk_name='chk_place_holder', key_words=''):
        '''
        generate gaussian 16 ONIOM route section. Base on settings in the config module.
        -------
        work_type   : ONIOM calculation type (support: spe, ...)
        chk_name    : filename of chk (QMMM.chk by default / self.name + _QMMM.chk in the default workflow use.)
        key_words   : allow additional key words
        -------
        support edit keywords directly in module Config
        '''
        chk = r'%chk='+chk_name+'.chk' + line_feed
        proc = '%nprocshared='+str(Config.n_cores) + line_feed
        mem = '%mem=' + str(Config.n_cores * Config.max_core) + 'MB' + line_feed
        if type(key_words) == str and key_words != '':
            keyword_line = '# '+' '.join(Config.Gaussian.keywords[work_type]+[key_words,]) + line_feed
        if type(key_words) == list:
            keyword_line = '# '+' '.join(Config.Gaussian.keywords[work_type]+key_words) + line_feed

        route = chk + proc + mem + keyword_line
        return route


    def _get_oniom_chrgspin(self, prmtop_path=None, spin_list=[1,1]):
        '''
        Determing charge and spin for each ONIOM layers. Base on *prmtop file* and layer settings in the *config* module.
        '''
        chrgspin = None
        # san check
        if prmtop_path == None:
            if self.prmtop_path == None:
                raise Exception('Please provide or use PDB2FF() to generate a prmtop file before PDB2QMMM')
            prmtop_path=self.prmtop_path

        # get charge list
        self.chrg_list_all = PDB.get_charge_list(prmtop_path)
        # init
        self.layer_chrgspin=[]
        for j in range(len(self.layer)):
            self.layer_chrgspin.append(float(0))
        # add charge to layers
        for i, chrg in enumerate(self.chrg_list_all):
            for j, layer in enumerate(self.layer):
                if i+1 in layer:
                    self.layer_chrgspin[j] += chrg

        # add spin
        if len(self.layer_chrgspin) != len(spin_list):
            raise Exception('spin specification need to match the layer setting. e.g.: spin_list=[h_spin, l_spin]')
        for i,spin in enumerate(spin_list):
            self.layer_chrgspin[i] = (self.layer_chrgspin[i], spin)    
        # make string
        if len(self.layer) == 2:
            c1 = str(round(self.layer_chrgspin[0][0]))
            s1 = str(round(self.layer_chrgspin[0][1]))
            c2 = str(round(self.layer_chrgspin[1][0]))
            s2 = str(round(self.layer_chrgspin[1][1]))
            c3 = c2
            s3 = s2
            chrgspin = ' '.join([c1,s1,c2,s2,c3,s3])
        else:
            raise Exception('Only support 2 layers writing charge and spin. Update in the future')

        if Config.debug >= 1:
            print(chrgspin)

        return chrgspin


    def _get_oniom_g16_coord(self):
        '''
        generate coordinate line. Base on *structure* and layer settings in the *config* module.
        Use element name as atom type for ligand atoms since they are mostly in QM regions.
        ---------------
        for a coord line:
            - element name
        	- atom name (from .lib)
        	- atom charge (from self.charge_list_all)
        	- freeze part (general option / some presupposition) 
        	- xyz (from self.stru)
        	- layer (general option / some presupposition)
        '''
        coord=''

        # amber default hold the chain - ligand - metal - solvent order
        a_id = 0
        for chain in self.stru.chains:
            for res in chain:
                for atom in res:
                    a_id += 1
                    # san check
                    if atom.id != a_id:
                        raise Exception('atom id error.')
                    if atom.id in self.layer[0]:
                        coord += atom.build_oniom('h', self.chrg_list_all[atom.id-1])
                    else:
                        # consider connection
                        cnt_info = None
                        repeat_flag = 0
                        for cnt_atom in atom.connect:
                            if cnt_atom.id in self.layer[0]:
                                if repeat_flag:
                                    raise Exception('A low layer atom is connecting 2 higher layer atoms')
                                cnt_info = ['H', cnt_atom.get_pseudo_H_type(atom), cnt_atom.id] 
                                repeat_flag = 1
                        # general low layer
                        coord += atom.build_oniom('l', self.chrg_list_all[atom.id-1], cnt_info=cnt_info)
        for lig in self.stru.ligands:
            for atom in lig:
                a_id += 1
                if atom.id != a_id:
                    raise Exception('atom id error.')
                if atom.id in self.layer[0]:
                    coord += atom.build_oniom('h', self.chrg_list_all[atom.id-1], if_lig=1)
                else:
                    if Config.debug >= 1:
                        print('\033[1;31;0m In PDB2QMMM in _get_oniom_g16_coord: WARNING: Found ligand atom in low layer \033[0m')
                    # consider connection
                    cnt_info = None
                    repeat_flag = 0 
                    for cnt_atom in atom.connect:
                        if repeat_flag:
                            raise Exception('A low layer atom is connecting 2 higher layer atoms')
                        if cnt_atom.id in self.layer[0]:
                            if Config.debug >= 1:
                                print('\033[1;31;0m In PDB2QMMM in _get_oniom_g16_coord: WARNING: Found ligand atom'+str(atom.id)+' in seperate layers \033[0m')
                            cnt_info = ['H', cnt_atom.get_pseudo_H_type(atom), cnt_atom.id]
                            repeat_flag = 1
                    coord += atom.build_oniom('l', self.chrg_list_all[atom.id-1], cnt_info=cnt_info, if_lig=1)
        for atom in self.stru.metalatoms:
            a_id += 1
            if atom.id != a_id:
                raise Exception('atom id error.')
            if atom.id in self.layer[0]:
                coord += atom.build_oniom('h', self.chrg_list_all[atom.id-1])
            else:
                coord += atom.build_oniom('l', self.chrg_list_all[atom.id-1])
        for sol in self.stru.solvents:
            for atom in sol:
                a_id += 1
                if atom.id != a_id:
                    raise Exception('atom id error.')
                if atom.id in self.layer[0]:
                    coord += atom.build_oniom('h', self.chrg_list_all[atom.id-1], if_sol=1)
                else:
                    # consider connection
                    cnt_info = None # for future update
                    repeat_flag = 0
                    for cnt_atom in atom.connect:
                        if cnt_atom.id in self.layer[0]:
                            if Config.debug >= 1:
                                print('\033[1;31;0m In PDB2QMMM in _get_oniom_g16_coord: WARNING: Found solvent atom'+str(atom.id)+' in seperate layers \033[0m')
                            if repeat_flag:
                                raise Exception('A low layer atom is connecting 2 higher layer atoms')
                            cnt_info = ['H', cnt_atom.get_pseudo_H_type(atom), cnt_atom.id] 
                            repeat_flag = 1
                    coord += atom.build_oniom('l', self.chrg_list_all[atom.id-1], cnt_info=cnt_info, if_sol=1)

        return coord


    def _get_oniom_g16_add_prm(self):
        '''
        Add missing parameters for protein and custom atom types
        1. addition parameters for metal element that not exist in ff96
        2. Commonly missing line for no reason: 'HrmBnd1    N   CT   HC     35.0000     109.5000'
        3. What if ligand or artificial residue appears in the low layer TODO
        4. missing parameters brought by the pseudo boundary H
        '''
        #1
        add_prm='HrmBnd1    N   CT   HC     35.0000     109.5000'+line_feed
        #2
        atom_rec=[]
        for atom in self.stru.metalatoms:
            if atom.parm != None:
                if atom.parm[0] not in atom_rec:
                    add_prm += 'VDW   '+ '   '.join(atom.parm)+line_feed
                    atom_rec.append(atom.parm[0])
        #3
        #TODO
        #4

        return add_prm

    @classmethod
    def get_charge_list(cls, prmtop_path):
        '''
        Get charge from the .prmtop file
        Take the (path) of .prmtop file and return a [list of charges] with corresponding to the atom sequence
        -----------------
        * Unit transfer in prmtop: http://ambermd.org/Questions/units.html
        '''
        with open(prmtop_path) as f:
            
            charge_list=[]
            line_index=0

            for line in f:
                
                line_index=line_index+1 #current line
                
                if line.strip() == r'%FLAG POINTERS':
                    format_flag=line_index
                if line.strip() == r'%FLAG CHARGE':
                    charge_flag=line_index

                if 'format_flag' in dir():
                    if line_index == format_flag+2:
                        N_atom=int(line.split()[0])
                        del format_flag
                        
                if 'charge_flag' in dir():
                    if line_index >= charge_flag+2 and line_index <= charge_flag+1+ceil(N_atom/5):
                        for i in line.strip().split():
                            charge_list.append(float(i)/18.2223)
        return charge_list

    '''
    ========
    MD Analysis 
    ========
    '''
    def nc2mdcrd(self, o_path='', point=None, start=1, end=-1, step=1, engine='cpptraj'):
        '''
        convert self.nc to a mdcrd file to read and operate.(self.nc[:-2]+'.mdcrd' by default)
        a easier way is to use pytraj directly.
        ---------------
        o_path: user assigned out path (self.nc[:-2]+'mdcrd' by default)
        point:  sample point. use value from self.conf_prod['nstlim'] and self.conf_prod['ntwx'] to determine step size.
        start:  start point
        end:    end point
        step:   step size
        engine: pytraj or cpptraj (some package conflict may cause pytraj not available)
        '''
        if self.nc == None:
            raise Exception('No nc file found. Please assign self.nc or run PDBMD first')
        else:
            if o_path == '':
                o_path= self.nc[:-2]+'mdcrd'
            if end == -1:
                end = 'last'
            if point != None:
                all_p = int(self.conf_prod['nstlim'])/int(self.conf_prod['ntwx'])
                step = int(all_p/point)

            if engine not in ['pytraj', 'cpptraj']:
                raise Exception('engine: pytraj or cpptraj')
            
            if engine == 'pytraj':
                pass

            if engine == 'cpptraj':
                cpp_in_path = self.cache_path+'/cpptraj_nc2mdcrd.in'
                cpp_out_path = self.cache_path+'/cpptraj_nc2mdcrd.out'
                with open(cpp_in_path,'w') as of:
                    of.write('parm '+self.prmtop_path+line_feed)
                    of.write('trajin '+self.nc+' '+str(start)+' '+end+' '+str(step)+line_feed)
                    of.write(f'autoimage{line_feed}') # this is important
                    of.write('trajout '+o_path+line_feed)
                    of.write('run'+line_feed)
                    of.write('quit'+line_feed)
                os.system('cpptraj -i '+cpp_in_path+' > '+cpp_out_path)

        self.mdcrd=o_path
        return o_path
            

    '''
    ========
    QM Cluster
    ========    
    '''
    def PDB2QMCluster(
        self, 
        atom_mask: str, 
        spin: int = 1, 
        o_dir: Union[str, None] = None, 
        tag: str = '', 
        QM: str = 'g16',
        g_route: Union[str, None] =None,
        ifchk: bool = 0, 
        val_fix: str = 'internal',
        if_cluster_job: bool = 1, 
        cluster: ClusterInterface = None,
        job_array_size: int = 0,
        period: int = 600,
        res_setting: Union[dict, str, None] = None,
        cpu_cores: Union[int, str, None] = None,
        cpu_mem: Union[int, str, None] = None,
        cluster_debug: bool = 0
    ) -> list :
        '''
        Build & Run QM cluster input from self.mdcrd with selected atoms according to atom_mask
        ---------
        Args:
        self (Pre-requisites):
            dir 
                - for *default dir*
            path 
                - for *get_stru()* 
                - need the same pdb representing the structure in the MD frames
            prmtop_path 
                - for get *charge* list 
                - need the prmtop file used as MD input
            mdcrd 
                - for *coordinates* of each QM cluster 
                - the mdcrd file that sampled from the traj
            prepi_path (val_fix='internal')
                - for get *connectivity (ligand part)* and fix free valances if they exist
                - the dict for prepin files for all ligands {'3_letter_name':'path_to_prepin_file', ...}
            (* the later 3 should of the consistancy from the same MD*)
        atom_mask:
            Amber-style atom mask for specifing the the QM region.
        spin: 
            specific spin state for the qm cluster. (default: 1)
        o_dir:
            The output directory that contain all QM input & output jobs.
        tag:
            Add tag to default out_dir. It will be {self.dir}/QM_cluster{tag}
        QM: 
            QM engine (default: g16)
        g_route: 
            Gaussian route line
        ifchk: 
            if save chk file of a gaussian job. (default: 0)
        val_fix: 
            fix free valance of truncated strutures
            = internal: add H to where the original connecting atom is. 
                        special fix for classical case:
                        "sele by residue" (cut N-C) -- adjust dihedral for added H on N.
            = openbabel TODO
        if_cluster_job:
            if submit the calculation to a cluster using the job manager
        cluster, job_array_size, period
            see Run_QM
        res_setting
            see Run_QM for detail (default: Config.Gaussian.QMCLUSTER_CPU_RES)
            provide specific keys you want to change the rest will remain default.
            (e.g. res_setting = {'node_cores':'8'} will still use a complete dictionary 
            with only node_cores modified.)
        cpu_cores, cpu_mem
            provide these for configure the gaussian input file when using a non-dict res_setting.
            these values should be the same as indicated in the res_setting.
        cluster_debug:
           1:  add also the qm cluster job obj to the pdb obj
        ---data---
        Attribute:
            self.frames
            self.qm_cluster_map (PDB atom id -> QM atom id)
        Return:
            self.qm_cluster_out
                the list of paths of the qm cluster output files
        '''
        # make folder
        if o_dir == None:
            o_dir = self.dir+'/QM_cluster'+tag
        mkdir(o_dir)
        # update stru
        self.get_stru()
        # get sele
        if val_fix == 'internal':
            sele_lines, sele_map = self.stru.get_sele_list(atom_mask, fix_end='H', prepi_path=self.prepi_path)
        else:
            sele_lines, sele_map = self.stru.get_sele_list(atom_mask, fix_end=None)
        self.qm_cluster_map = sele_map
        # get chrgspin
        chrgspin = self._get_qmcluster_chrgspin(sele_lines, spin=spin)
        if Config.debug >= 1:
            print('Charge: '+str(chrgspin[0])+' Spin: '+str(chrgspin[1]))
        # get res setting
        # overwrite those if cluster job is true
        if if_cluster_job:
            # default values TODO should contain them in Config.Gaussian
            res_setting = type(self)._get_default_res_setting_qmcluster(res_setting)
            if QM in ['g16','g09']:
                if isinstance(res_setting, dict):
                    cpu_cores = res_setting['node_cores']
                    cpu_mem = res_setting['mem_per_core']
                elif cpu_cores is None or cpu_mem is None:
                    raise TypeError('ERROR: Requires cpu_cores and cpu_mem input for configuring Gaussian input file if using CPU and provide res_setting not as a dict')
                cpu_mem = int(float(cpu_mem.rstrip('GB')) * 1024) # convert to number in MB
        else:
            cpu_cores = Config.n_cores
            cpu_mem = Config.max_core

        #make inp files
        frames = Frame.fromMDCrd(self.mdcrd)
        self.frames = frames
        if QM in ['g16','g09']:
            gjf_paths = []
            if Config.debug >= 1:
                print('Writing QMcluster gjfs.')
            for i, frame in enumerate(frames):
                gjf_path = o_dir+'/qm_cluster_'+str(i)+'.gjf'
                frame.write_sele_lines(sele_lines, out_path=gjf_path, g_route=g_route, g_cores=cpu_cores, g_mem_cores=cpu_mem, chrgspin=chrgspin, ifchk=ifchk)
                gjf_paths.append(gjf_path)
            # Run inp files
            if if_cluster_job:
                Run_QM_out = PDB.Run_QM( gjf_paths, 
                                         prog=QM, 
                                         if_cluster_job=if_cluster_job, 
                                         cluster = cluster,
                                         job_array_size = job_array_size,
                                         period = period,
                                         res_setting=res_setting,
                                         cluster_debug=cluster_debug)
                if cluster_debug:
                    qm_cluster_out_paths = Run_QM_out[0]
                    self.qm_cluster_jobs = Run_QM_out[1]
                else:
                    qm_cluster_out_paths = Run_QM_out
            else:
                qm_cluster_out_paths = PDB.Run_QM(gjf_paths, prog=QM, if_cluster_job=0)
            # get chk files if ifchk
            if ifchk:
                qm_cluster_chk_paths = []
                for gjf in gjf_paths:
                    chk_path = gjf[:-3]+'chk'
                    qm_cluster_chk_paths.append(chk_path)
        if QM=='ORCA':
            pass

        self.qm_cluster_out = qm_cluster_out_paths
        if ifchk:
            self.qm_cluster_chk = qm_cluster_chk_paths
            return self.qm_cluster_out, self.qm_cluster_chk
        
        return self.qm_cluster_out
    

    def _get_qmcluster_chrgspin(self, sele, spin=1):
        '''
        get charge for qmcluster of sele from self.prmtop
        # TODO refine charge count for qm region cut interface
        '''
        # get chrg list
        chrg_list_all = PDB.get_charge_list(self.prmtop_path)
        # sum with sele
        sele_chrg = 0
        for sele_atom in sele.keys():
            if '-' in sele_atom:
                continue
            # skip backbone atoms (cut CA-CB when calculate charge)
            if sele_atom[-1] == 'b':
                continue

            # clean up
            if sele_atom[-1] not in '1234567890':
                sele_id = sele_atom[:-1]
            else:
                sele_id = sele_atom
            sele_chrg += chrg_list_all[int(sele_id)-1]

        return (round(sele_chrg), spin)

    @staticmethod
    def _get_default_res_setting_qmcluster(res_setting):
        if res_setting is None:
            res_setting = dict()
        if isinstance(res_setting, dict):
            res_setting_holder = copy.deepcopy(Config.Gaussian.QMCLUSTER_CPU_RES)
            # replace assigned key in default dict
            for k, v in res_setting.items():
                res_setting_holder[k] = v
            return res_setting_holder
        if isinstance(res_setting, str):
            return res_setting

    @classmethod
    def Run_QM(
        cls, 
        inp: list[str], 
        prog: str = 'g16', 
        if_cluster_job: bool = 1,
        cluster: ClusterInterface = None,
        job_array_size: int = 0,
        period: int = 600,
        res_setting: dict = None,
        clean_job_cluster_log: bool = True,
        cluster_debug: bool = 0
    ):
        '''
        Run QM with {prog} for {inp} files and return paths of output files.
        Args:
        inp: 
            list of paths of input files
        prog:
            the QM program (default: g16)
        if_cluster_job:
            if submit the QM calculation to a HPC (default: 1)
        cluster:
            The cluster used when if_cluster_job is 1.
        job_array_size:
            how many jobs are allowed to submit simultaneously. (default: 0 means len(inp))
            (e.g. 5 for 100 jobs means run 20 groups. All groups will be submitted and 
            in each group, submit the next job only after the previous one finishes.)
        period:
            the time cycle for each job state change (Unit: s) (default: 600)
        res_setting:
            resource setting dictionary
        cluster_debug:
            1: return also the job objects
            0: return only the out file paths

        TODO put this individually as part of the qm interface
             maybe introduct the current executor object to decouple this module with the job manager.
        '''
        if if_cluster_job:
            #san check
            if not isinstance(cluster, ClusterInterface):
                raise TypeError('cluster job need a cluster (ClusterInterface object) input')
            
            if prog == 'g16':
                outs = []
                # config jobs
                jobs = []
                for gjf_path in inp:
                    out_path = gjf_path.removesuffix('gjf')+'out'
                    jobs.append(cls._make_single_g16_job(gjf_path, out_path, cluster, res_setting))
                    outs.append(out_path)
                # submit and run in array
                if Config.debug > 0:
                    print(f'''Running QM array on {cluster.NAME}: number: {len(jobs)} size: {job_array_size} period: {period}''')
                job_manager.ClusterJob.wait_to_array_end(jobs, period, job_array_size)
                if clean_job_cluster_log:
                    target_dir = f"{os.path.dirname(inp[0])}/cluster_log/"          
                    mkdir(target_dir)
                    for job in jobs:
                        job: job_manager.ClusterJob
                        new_path = shutil.move(job.job_cluster_log, target_dir)
                        job.job_cluster_log = new_path
                if cluster_debug:
                    return outs, jobs
                return outs
        else:
            # local job
            if prog == 'g16':
                outs = []
                for gjf in inp:
                    out = gjf[:-3]+'out'
                    if Config.debug > 1:
                        print('running: '+Config.Gaussian.g16_exe+' < '+gjf+' > '+out)
                    os.system(Config.Gaussian.g16_exe+' < '+gjf+' > '+out)
                    outs.append(out)
                return outs

            if prog == 'g09':
                outs = []
                for gjf in inp:
                    out = gjf[:-3]+'out'
                    if Config.debug > 1:
                        print('running: '+Config.Gaussian.g09_exe+' < '+gjf+' > '+out)
                    os.system(Config.Gaussian.g09_exe+' < '+gjf+' > '+out)
                    outs.append(out)
                return outs

    @classmethod
    def _make_single_g16_job(
            cls, 
            gjf_path: str, 
            out_path: str, 
            cluster: ClusterInterface,
            res_setting: dict
        ) -> job_manager.ClusterJob :
        '''
        job for submit g16 for gjf > out to cluster
        return a ClusterJob object
        '''
        cmd = f'{Config.Gaussian.g16_exe} < {gjf_path} > {out_path}'
        # interface check
        if 'G16_ENV' not in dir(cluster):
            raise Exception('RunQM(prog = g16) requires the input cluster have the G16_ENV attr')
        job = job_manager.ClusterJob.config_job(
            commands = cmd,
            cluster = cluster,
            env_settings = cluster.G16_ENV['CPU'],
            res_keywords = res_setting,
            sub_dir = './', # because gjf path are relative
            sub_script_path = gjf_path.removesuffix('gjf')+'cmd'
        )
        return job

    def get_fchk(self, keep_chk=0):
        '''
        transfer Gaussian chk files to fchk files using formchk
        ----------
        keep_chk: if not delete original chk file (default: 0)
        '''
        # san check
        if len(self.qm_cluster_chk) == 0:
            raise Exception('No chk file in self.qm_cluster_chk.')

        # formchk
        fchk_paths = []
        for chk in self.qm_cluster_chk:
            fchk = chk[:-3]+'fchk'
            if Config.debug > 1:
                print('running: '+'formchk '+chk+' '+fchk)
            run('formchk '+chk+' '+fchk, check=True, text=True, shell=True, capture_output=True)
            fchk_paths.append(fchk)
            # keep chk
            if not keep_chk:
                if Config.debug > 1:
                    print('removing: '+chk)
                os.remove(chk)

        self.qm_cluster_fchk = fchk_paths
        return self.qm_cluster_fchk
        
    def gaussian_error_handling():
        '''
        TODO handle gaussian errors. can be enabled in RunQM
        use maps in Config.Gaussian
        '''
        pass
    '''
    ========
    Analysis 
    ========
    '''
    def get_field_strength(self, atom_mask, a1=None, a2=None, bond_p1='center', p1=None, p2=None, d1=None):
        '''
        use frame coordinate from *mdcrd* and MM charge from *prmtop* to calculate the field strength of *p1* along *p2-p1* or *d1*
        atoms in *atom_mask* is included. (TODO: or an exclude one?)
        -------------------------------------
        a1 a2:  id of atoms compose the bond
        bond_p1:method to generate p1
                - center
                - a1
                - TODO exact dipole center of the current structure
                - ...
        p1:     the point where E is calculated
        p2:     a point to fix d1
        d1:     the direction E is projected
        return an ensemble field strengths
        '''
        Es = []
        chrg_list = PDB.get_charge_list(self.prmtop_path)

        # san check
        if a1 == None and p1 == None:
            raise Exception('Please provide a 1nd atom (a1=...) or point (p1=...) where E is calculated ')
        if a2 == None and p2 == None and d1 == None:
            raise Exception('Please provide a 2nd atom (a2=...) or point (p2=...) or a direction (d1=...)')
        if a1 == None and a2 == None and bond_p1 == 'center':
            raise Exception('Please provide a both atom (a1=... a2=...) when bond_p1 = center; Or change bond_p1 to a1 to calculate E at a1')
        if a1 != None and a2 != None and bond_p1 not in ['center', 'a1']:
            raise Exception('Only support p1 selection in center or a1 now')

        if self.frames == None:
            self.frames = Frame.fromMDCrd(self.mdcrd)

        # decode atom mask (stru corresponding to mdcrd structures)
        atom_list = decode_atom_mask(self.stru, atom_mask)

        for frame in self.frames:
            #get p2
            if a2 != None:
                p2 = frame.coord[a2-1]
            #get p1
            if a1 != None:
                if bond_p1 == 'a1':
                    p1 = frame.coord[a1-1]
                if bond_p1 == 'center':
                    p1 = get_center(frame.coord[a1-1], p2)
                if bond_p1 == 'xxx':
                    pass
            # sum up field strength
            E = 0
            for atom_id in atom_list:
                # search for coord and chrg
                coord = frame.coord[atom_id-1]
                chrg = chrg_list[atom_id-1]
                E += get_field_strength_value(coord, chrg, p1, p2=p2, d1=d1)
            Es.append(E)

        return Es

    @classmethod
    def get_bond_dipole(cls, qm_fch_paths, a1, a2, prog='Multiwfn'):
        '''
        get bond dipole using wfn analysis with fchk files.
        -----------
        Args:
            qm_fch_paths: paths of fchk files 
                        * requires correponding out files with only ext difference
                        * (if want to compare resulting coord to original mdcrd/gjf stru)
                            requires nosymm in gaussian input that generate the fch file.
            a1          : QM I/O id of atom 1 of the target bond
            a2          : QM I/O id of atom 2 of the target bond
            prog        : program for wfn analysis (default: multiwfn)
                          **Multiwfn workflow**
                            1. Multiwfn xxx.fchk < parameter_file > output
                            the result will be in ./LMOdip.txt 
                            2. extract value and project to the bond accordingly
        Returns:
            Dipoles     : A list of dipole data in a form of [(dipole_norm_signed, dipole_vec), ...]
                          *dipole_norm_signed* is the signed norm of the dipole according to its projection
                                               to be bond vector.
                          *dipole_vec* is the vector of the dipole
        -----------
        LMO bond dipole: 
        (Method: Multiwfn manual 3.22/4.19.4)
        2-center LMO dipole is defined by the deviation of the eletronic mass center relative to the bond center.
        Dipole positive Direction: negative(-) to positive(+).
        Result direction: a1 -> a2
        
        REF: Lu, T.; Chen, F., Multiwfn: A multifunctional wavefunction analyzer. J. Comput. Chem. 2012, 33 (5), 580-592.
        ''' # TODO support array job
        Dipoles = []

        if prog == 'Multiwfn':
            # self.init_Multiwfn()
            ref_path = qm_fch_paths[0]
            mltwfn_in_path = ref_path[:-(len(ref_path.split('.')[-1])+1)]+'_dipole.in'

            with open(mltwfn_in_path, 'w') as of:
                of.write('19'+line_feed)
                of.write('-8'+line_feed)
                of.write('1'+line_feed)
                of.write('y'+line_feed)
                of.write('q'+line_feed)

            bond_id_pattern = r'\( *([0-9]+)[A-Z][A-z]? *- *([0-9]+)[A-Z][A-z]? *\)'
            bond_data_pattern = r'X\/Y\/Z: *([0-9\.\-]+) *([0-9\.\-]+) *([0-9\.\-]+) *Norm: *([0-9\.]+)'
            
            for fchk in qm_fch_paths:
                # get a1->a2 vector from .out (update to using fchk TODO)
                G_out_path = fchk[:-len(fchk.split('.')[-1])]+'out'
                with open(G_out_path) as f0:
                    coord_flag = 0
                    skip_flag = 0
                    for line0 in f0:
                        if 'Input orientation' in line0:
                            coord_flag = 1
                            continue
                        if coord_flag:
                            # skip first 4 lines
                            if skip_flag <= 3:
                                skip_flag += 1
                                continue
                            if '-------------------------' in line0:
                                break
                            l_p = line0.strip().split()
                            if str(a1) == l_p[0]:
                                coord_a1 = np.array((float(l_p[3]), float(l_p[4]), float(l_p[5])))
                            if str(a2) == l_p[0]:
                                coord_a2 = np.array((float(l_p[3]), float(l_p[4]), float(l_p[5])))
                Bond_vec = (coord_a2 - coord_a1)
                
                # Run Multiwfn
                mltwfn_out_path = fchk[:-len(fchk.split('.')[-1])]+'dip'
                if Config.debug >= 2:
                    print('Running: '+Config.Multiwfn.exe+' '+fchk+' < '+mltwfn_in_path)
                run(Config.Multiwfn.exe+' '+fchk+' < '+mltwfn_in_path, check=True, text=True, shell=True, capture_output=True)
                run('mv LMOdip.txt '+mltwfn_out_path, check=True, text=True, shell=True, capture_output=True)
                run('rm LMOcen.txt new.fch', check=True, text=True, shell=True, capture_output=True)
                
                # get dipole
                with open(mltwfn_out_path) as f:
                    read_flag = 0
                    for line in f:
                        if line.strip() == 'Two-center bond dipole moments (a.u.):':
                            read_flag = 1
                            continue
                        if read_flag:
                            if 'Sum' in line:
                                raise Exception('Cannot find bond:'+str(a1)+'-'+str(a2)+line_feed)
                            Bond_id = re.search(bond_id_pattern, line).groups()
                            # find target bond
                            if str(a1) in Bond_id and str(a2) in Bond_id:
                                Bond_data = re.search(bond_data_pattern, line).groups()
                                dipole_vec = (float(Bond_data[0]) ,float(Bond_data[1]) ,float(Bond_data[2]))
                                # determine sign
                                if np.dot(np.array(dipole_vec), Bond_vec) > 0:
                                    dipole_norm_signed = float(Bond_data[3])
                                else:
                                    dipole_norm_signed = -float(Bond_data[3])
                                break
                Dipoles.append((dipole_norm_signed, dipole_vec))

        return Dipoles

    @classmethod
    def init_Multiwfn(cls, n_cores=None):
        '''
        initiate Multiwfn with settings in Config
        '''
        # set nthreads
        if n_cores == None:
            n_cores = str(Config.n_cores)
        if Config.debug >= 1:
            print("Running: "+"sed -i 's/nthreads= *[0-9][0-9]*/nthreads=  "+n_cores+"/' "+Config.Multiwfn.DIR+"/settings.ini")
        run("sed -i 's/nthreads= *[0-9][0-9]*/nthreads=  "+n_cores+"/' "+Config.Multiwfn.DIR+"/settings.ini", check=True, text=True, shell=True, capture_output=True)

    def get_rosetta_ddg(self, rosetta_home: str, muta_groups: List[tuple], relaxed_pdb: str,
                        niter: int = 10,
                        ddg_dir: str = None,
                        in_place: bool = False, 
                        if_cluster_job: bool = True,
                        cluster: ClusterInterface = None,
                        period: int = 120,
                        job_array_size: int = 100,
                        res_setting: Union[dict, None] = None,
                        cluster_debug: bool = 0 ) -> dict:
        '''
        MVP function to obtain Rosetta ddG stability score for current structure in PDB()
        ref: https://www.rosettacommons.org/docs/latest/cartesian-ddG
        due to the fact that cartesian_ddg can only benefit from 1 cpu core, we will parallel
        multi-mutants with ARMer job array to speed up.
        '''
        TEMP_ROSETTA_RES = {    'core_type' : 'cpu',
                                'nodes':'1',
                                'node_cores' : '1',
                                'job_name' : 'EnzyHTP_ddG',
                                'partition' : 'production',
                                'mem_per_core' : '2G',
                                'walltime' : '1-00:00:00',
                                'account' : 'yang_lab_csb'}
        if ddg_dir is None:
            ddg_dir = f'{self.dir}/ddg/'
        mkdir(ddg_dir)
        ddg_jobs = []
        ddg_result_files = {}

        for i, mut_group in enumerate(muta_groups):
            group_dir = f'{ddg_dir}group_{i}/'
            mkdir(group_dir)
            flag_file_path = f'{group_dir}/ddg_options.txt'
            muta_file_path = f'{group_dir}/mutation.txt'
            # make option file
            with open(flag_file_path, 'w') as of:
                of.write(f'''-in:file:s {os.path.abspath(relaxed_pdb)}
-ddg::mut_file ./{os.path.relpath(muta_file_path, group_dir)}
-ddg:iterations {niter}
-force_iterations false
-ddg::score_cutoff 1.0
-ddg::cartesian
-ddg::dump_pdbs false
-fa_max_dis 9.0
-score:weights talaris2014_cart
-ddg::legacy false
-mute all''')
            # make mutant file
            with open(muta_file_path, 'w') as of: # make this into different jobs (different folder containing options mutations submission)
                of.write(f'total {len(mut_group)}{line_feed}')
                of.write(f'{len(mut_group)}{line_feed}')
                for r_mutaflag in mut_group:
                    of.write(f'{r_mutaflag}{line_feed}')
            # make job
            ddg_cmd = f'{rosetta_home}/source/bin/cartesian_ddg.static.linuxgccrelease @ ./{os.path.relpath(flag_file_path, group_dir)}'
            ddg_job = job_manager.ClusterJob.config_job(
                        commands = ddg_cmd,
                        cluster = cluster,
                        env_settings = '',
                        res_keywords = TEMP_ROSETTA_RES,
                        sub_dir = group_dir,
                        sub_script_path = f'{group_dir}/submit_ddg.cmd')
            ddg_jobs.append(ddg_job)
            ddg_result_files[mut_group] = f'{group_dir}/mutation.ddg'
        
        job_manager.ClusterJob.wait_to_array_end(ddg_jobs, period, job_array_size)
        # support multi groups?
        result = {}
        for mutants, ddg_file in ddg_result_files.items():
            ddg_value = type(self).get_rosetta_ddg_result(ddg_file, niter)
            result[mutants] = ddg_value
        return result
    
    def relax_with_rosetta(self, rosetta_home: str, nstruct_relax: int = 20,
                            in_place: bool = False, if_clean: bool = True,
                            if_cluster_job: bool = True,
                            cluster: ClusterInterface = None,
                            period: int = 120,
                            res_setting: Union[dict, None] = None,
                            cluster_debug: bool = 0 ) -> str:
        '''
        MVP function to relax the structure with Rosetta
        only can work with ddG calculation right now
        ref: https://www.rosettacommons.org/docs/latest/application_documentation/structure_prediction/relax
        '''
        TEMP_ROSETTA_RES = {    'core_type' : 'cpu',
                                'nodes':'1',
                                'node_cores' : '24',
                                'job_name' : 'EnzyHTP_r_relax',
                                'partition' : 'production',
                                'mem_per_core' : '2G',
                                'walltime' : '1-00:00:00',
                                'account' : 'yang_lab_csb'}

        int_pdb_path_1 = f'{self.path.removesuffix(".pdb")}_rosetta.pdb'
        int_pdb_path_2 = f'{self.path.removesuffix(".pdb")}_relaxed.pdb'
        # remove ligands & clean for rosetta (tmp)
        run(f'{rosetta_home}/tools/protein_tools/scripts/clean_pdb.py --allchains {self.path} ignorechain', # TODO replace this as described in https://www.shaoqz.cn/2022/02/05/Rosetta-Notes/
            check=True, text=True, shell=True, capture_output=True)
        run(f'mv {self.path.split("/")[-1].removesuffix(".pdb")}_ignorechain.pdb {int_pdb_path_1}',
            check=True, text=True, shell=True, capture_output=True)        
        run(f'rm {self.path.split("/")[-1].removesuffix(".pdb")}_ignorechain.fasta',
            check=True, text=True, shell=True, capture_output=True)        
        # relax
        relax_cmd = f'mpiexec -np {TEMP_ROSETTA_RES["node_cores"]} {rosetta_home}/source/bin/relax.mpi.linuxgccrelease -s {int_pdb_path_1} -use_input_sc -ignore_unrecognized_res -nstruct {nstruct_relax} -fa_max_dis 9.0'
        relax_job = job_manager.ClusterJob.config_job(
                    commands = relax_cmd,
                    cluster = cluster,
                    env_settings = ['module load GCC/5.4.0-2.26', 'module load OpenMPI'],
                    res_keywords = TEMP_ROSETTA_RES,
                    sub_dir = self.dir,
                    sub_script_path = f'{self.dir}/submit_relax.cmd')
        relax_job.submit()
        relax_job.wait_to_end(period=period)
        # get relaxed pdb
        # get the lowest score
        target_idx = type(self).get_rosetta_lowest_score(self.dir+'/score.sc')
        relaxed_pdb = f'{self.path.removesuffix(".pdb")}_rosetta_{target_idx:0>4}.pdb'
        run(f'cp {relaxed_pdb} {int_pdb_path_2}',
            check=True, text=True, shell=True, capture_output=True)
        
        if if_clean:
            for f_p in map(lambda x: f'{self.path.removesuffix(".pdb")}_rosetta_{x:0>4}.pdb', range(1,nstruct_relax+1)):
                os.remove(f_p)
            os.remove(self.dir+'/score.sc')
        return int_pdb_path_2
    
    @classmethod
    def get_rosetta_lowest_score(cls, score_file):
        '''
        get the structure index with the lowest score in the score.sc file
        '''
        score_df = pd.read_csv(score_file, delim_whitespace=True, header=1)
        min_sc_idx = score_df['total_score'].idxmin() + 1
        return min_sc_idx

    @classmethod
    def get_rosetta_ddg_result(cls, ddg_file: str, niter: int):
        '''
        get the result ddg from the .ddg file
        '''
        score_df = pd.read_csv(ddg_file, delim_whitespace=True, header=None)
        wt_sc_mean = score_df[3][:niter].mean()
        mut_sc_mean = score_df[3][niter:].mean()
        return mut_sc_mean - wt_sc_mean

    @classmethod
    def get_rmsd(cls, prmtop_path: str, traj_path: str, mask: str) -> float:
        """
        mvp function for RMSD calculation
        """
        tmp_rmsd_in = "./cpptraj_rmsd.in"
        tmp_rmsd_out = "./rmsd.dat"

        with open(tmp_rmsd_in, "w") as of:
            of.write(f"""parm {prmtop_path}
trajin {traj_path}
autoimage
rmsd {mask} first mass
average crdset AVE {mask}
run
autoimage
rmsd {mask} ref AVE * out {tmp_rmsd_out} mass
run
""")
        try:
            run(f"cpptraj -i {tmp_rmsd_in}",
            check=True, text=True, shell=True, capture_output=True)
        except CalledProcessError as e:
            print(e.stdout, e.stderr, sep=os.linesep)
            sys.exit(1)
        result = cls.get_cpptraj_rmsd_result(tmp_rmsd_out)
        os.remove(tmp_rmsd_in)
        os.remove(tmp_rmsd_out)
        return result

    @classmethod
    def get_cpptraj_rmsd_result(cls, result_path: str) -> float:
        """
        mvp function for RMSD calculation
        """
        result_df = pd.read_csv(result_path, delim_whitespace=True)
        return result_df.iloc[:, 1].mean()
    
    @classmethod
    def get_sasa_ratio(cls, 
                 prmtop_path: str, traj_path: str,
                 mask_pro: str, mask_pro_target: str, mask_sub: str,
                 tmp_dir: str = '.') -> float:
        """
        mvp function for SASA calculation
        Args:
            prmtop_path
            traj_path
            mask_pro: mask selection for the protein
            mask_pro_target: mask selection for the subsection in the protein as the SASA target
            mask_sub: mask selection for the substrate
            tmp_dir: temporary dir for divided traj files
        Return:
            (sasa_sub/sasa_pro).mean() average sasa ratio from each frame
        # WARNING: @biopaul919 reported mt.load have a bug that make the result wrong.
        """
        import mdtraj as mt
        # divide traj
        traj_parm_1, traj_parm_2 = cls.divide_traj(prmtop_path, traj_path, mask_pro, mask_sub, tmp_dir)
        # protein sasa
        traj_1 = mt.load(traj_parm_1[0], top=traj_parm_1[1])
        sasa_df_pro = mt.shrake_rupley(traj_1, n_sphere_points=5000, mode='residue')
        target_resi_idx = map(lambda x: int(x.strip()), mask_pro_target[1:].split(','))
        holder = np.zeros(len(sasa_df_pro)) # holder for sasa_pro of each frame
        test = 0
        for idx in target_resi_idx:
            test +=  sasa_df_pro[0][idx-1]
            holder += sasa_df_pro[:, idx-1]
        sasa_pro_by_frame = holder*100
        # substrate sasa
        traj_2 = mt.load(traj_parm_2[0], top=traj_parm_2[1])
        sasa_df_sub = mt.shrake_rupley(traj_2, n_sphere_points=5000, mode='residue')
        sasa_sub_by_frame = sasa_df_sub[:, 0]*100

        # clean up
        os.remove(traj_parm_1[0])
        os.remove(traj_parm_1[1])
        os.remove(traj_parm_2[0])
        os.remove(traj_parm_2[1])

        return (sasa_sub_by_frame/sasa_pro_by_frame).mean()

    @classmethod
    def get_surface_ratio(cls, 
                 prmtop_path: str, traj_path: str,
                 mask_pro: str, mask_pro_target: str, mask_sub: str, dynamic_sele: bool=False,
                 dot_solvent: int = 0, dot_density: int = 2,  
                 tmp_dir: str = '.', traj_start: str=None, traj_end: str=None,
                 if_complete_data: bool=False) -> float:
        """
        mvp function for SES calculation
        Args:
            prmtop_path
            traj_path
            mask_pro: mask selection for the protein
            mask_pro_target: mask selection for the subsection in the protein as the SASA target
            mask_sub: mask selection for the substrate
            dynamic_sele: if update selection every frame, must use if distance based pattern is used
                          in mask_pro_target
            dot_solvent: surface type (0: ses 1: sasa)
            dot_density: surface qulity
            tmp_dir: temporary dir for divided traj files
        Return:
            avg(ses_sub)/avg(ses_pro)
        """
        import pymol2
        if dynamic_sele:
            # 1. split states frist
            # 2. iterate using state names
            # 3. use pattern like f"(br. (((({obj_name} & {mask}) around 4) and {obj_name}) and (not n. C+H+O+N) and (not elem H)) and (not resi 1449+1181)) & n. CA"
            # ref https://github.com/shaoqx/CalcKit/tree/master/pymol_script/anchor_detect.py
            sys.exit(1) # TODO place holder
        # divide traj
        traj_parm_1, traj_parm_2 = cls.divide_traj(
            prmtop_path, traj_path, mask_pro, mask_sub, 
            tmp_dir, traj_start, traj_end)
        # protein ses
        ## start pymol
        p_session = pymol2.PyMOL()
        p_session.start()
        p_session.cmd.set("dot_solvent", dot_solvent)
        p_session.cmd.set("dot_density", dot_density)
        if Config.debug < 2:
            p_session.cmd.feedback("disable","all","everything")
        ## load traj protein
        pro_traj_obj = p_session.cmd.get_unused_name("pro_traj")
        p_session.cmd.load(traj_parm_1[1], pro_traj_obj)
        p_session.cmd.load_traj(traj_parm_1[0], pro_traj_obj)
        ## define selection
        if mask_pro_target.startswith(":"):
            target_resi_idx = map(lambda x: str(x.strip()), mask_pro_target[1:].split(','))
            pymol_sele_pattern = f"resi {'+'.join(target_resi_idx)}"
        else:
            pymol_sele_pattern = mask_pro_target
        ## calculate SES
        num_frames = p_session.cmd.count_states(pro_traj_obj)
        ses_pro_by_frame = []
        for i in range(num_frames):
            ses_pro_by_frame.append(
                p_session.cmd.get_area(selection=f"{pymol_sele_pattern} & {pro_traj_obj}", state=i+1)
            )
        p_session.cmd.delete(pro_traj_obj)
        # substrate sasa
        ## load traj substrate
        sub_traj_obj = p_session.cmd.get_unused_name("sub_traj")
        p_session.cmd.load(traj_parm_2[1], sub_traj_obj)
        p_session.cmd.load_traj(traj_parm_2[0], sub_traj_obj)
        ## calculate SES
        num_frames = p_session.cmd.count_states(sub_traj_obj)
        ses_sub_by_frame = []
        for i in range(num_frames):
            ses_sub_by_frame.append(
                p_session.cmd.get_area(selection=sub_traj_obj, state=i+1)
            )
        p_session.cmd.delete(sub_traj_obj)
        # clean up
        os.remove(traj_parm_1[0])
        os.remove(traj_parm_1[1])
        os.remove(traj_parm_2[0])
        os.remove(traj_parm_2[1])
        p_session.stop()
        if if_complete_data:
            return np.array(ses_sub_by_frame).mean(), np.array(ses_pro_by_frame).mean()
        return (np.array(ses_sub_by_frame)/np.array(ses_pro_by_frame)).mean()

    @classmethod
    def divide_traj(cls, prmtop_path: str, traj_path: str, mask1: str, mask2: str,
                    tmp_dir: str = '.', 
                    traj_start: str=None, traj_end: str=None) -> tuple:
        """
        divide traj into 2 trajs defined by mask1 and mask2
        """
        tmp_cpptraj_in = f'{tmp_dir}/cpptraj_divide.in'
        tmp_traj_1 = f'{tmp_dir}/traj_1.nc'
        tmp_traj_2 = f'{tmp_dir}/traj_2.nc'
        prmtop_name = prmtop_path.split('/')[-1]
        tmp_prmtop_1 = f'{tmp_dir}/traj_1.{prmtop_name}'
        tmp_prmtop_2 = f'{tmp_dir}/traj_2.{prmtop_name}'
        start_end_pattern = ''
        if traj_start and traj_end is not None:
            start_end_pattern = f' {traj_start} {traj_end}'

        with open(tmp_cpptraj_in, 'w') as of:
            of.write(f'''parm {prmtop_path}
trajin {traj_path}{start_end_pattern}
autoimage
strip !({mask1}) nobox outprefix {tmp_dir}/traj_1
outtraj {tmp_traj_1}
unstrip
strip !({mask2}) nobox outprefix {tmp_dir}/traj_2
outtraj {tmp_traj_2}
''')
        try:
            run(f"cpptraj -i {tmp_cpptraj_in}",
                check=True, text=True, shell=True, capture_output=True)
        except CalledProcessError as e:
            print(e.stdout, e.stderr, sep=os.linesep)
            sys.exit(1)
        os.remove(tmp_cpptraj_in)
        return (tmp_traj_1, tmp_prmtop_1), (tmp_traj_2, tmp_prmtop_2)

    def get_residue_pka(self, target_mask: str) -> float:
        '''
        mvp function for pka calculation for a specific residue using PROPKA
        Args:
            target_mask: pseudo AMBER style. e.g: ':1,2,3'
        '''
        target_resis = map(lambda x: int(x.strip()), 
            target_mask.removeprefix(':').split(','))
        # run propka (TODO go to interface)
        from propka.lib import loadOptions
        from propka.input import read_parameter_file, read_molecule_file
        from propka.parameters import Parameters
        from propka.molecular_container import MolecularContainer
        options = loadOptions([self.path]) # use default in mvp
        pdbfile = options.filenames[0]
        parameters = read_parameter_file(options.parameters, Parameters())
        my_molecule = MolecularContainer(parameters, options)
        my_molecule = read_molecule_file(pdbfile, my_molecule)
        my_molecule.calculate_pka()
        conformation = my_molecule.conformations['AVR']
        residues_pka = {}
        for group in conformation.groups:
            row_dict = {}
            atom = group.atom
            row_dict["res_num"] = atom.res_num
            row_dict["ins_code"] = atom.icode
            row_dict["res_name"] = atom.res_name
            row_dict["chain_id"] = atom.chain_id
            row_dict["group_label"] = group.label
            row_dict["group_type"] = getattr(group, "type", None)
            row_dict["pKa"] = group.pka_value
            row_dict["model_pKa"] = group.model_pka
            row_dict["buried"] = group.buried
            if group.coupled_titrating_group:
                row_dict["coupled_group"] = group.coupled_titrating_group.label
            else:
                row_dict["coupled_group"] = None
            residues_pka[atom.res_num] = row_dict
        return list(map(lambda y: residues_pka[y]['pKa'], target_resis))
    
    def get_mmpbsa_binding(
        self, 
        ligand_mask: str, 
        igb: int = 5,
        in_file: str = None,
        overwrite_in: str = 0,
        if_cluster_job: bool = True,
        cluster: ClusterInterface = None,
        period: int = 30,
        res_setting: Union[dict, None] = None,
        cluster_debug: bool = 0):
        """mvp function for getting the MMPBSA binding assessment for
        {ligand_mask} from the trajectory"""
        traj_file = self.mdcrd
        dr_prmtop, dl_prmtop, dc_prmtop, sc_prmtop = self.make_mmpbsa_prmtops(ligand_mask, igb)
        mmpbsa_out_file = self.run_mmpbsa(
            dr_prmtop, dl_prmtop, dc_prmtop, sc_prmtop, traj_file,
            in_file, overwrite_in, 
            if_cluster_job, cluster, period, res_setting, cluster_debug
            )
        result = type(self).extract_mmpbsa_out(mmpbsa_out_file)

        #clean
        if Config.debug < 2:
            os.remove(dr_prmtop)
            os.remove(dl_prmtop)
            os.remove(dc_prmtop)
            os.remove(sc_prmtop)
        if Config.debug < 1:    
            os.remove(mmpbsa_out_file)

        return result

    def make_mmpbsa_prmtops(self, ligand_mask: str, igb: int=5, use_ante_mmpbsa: bool=True):
        '''
        Make prmtop files for the MMPB(GB)SA calculation
        1. make new pdbs
        2. use tLeap to generate these from pdbs (use existing frcmod file)
        Args:
            frag_str:
                A amber like str to define two fragments
            igb:
                gb method used
        Returns:
            dr_prmtop, dl_prmtop, dc_prmtop, sc_prmtop 
            for dry receptor, dry ligand, dry complex, and solvate complex, respectively
        '''
        temp_dir = f"{self.dir}/temp/"
        mkdir(temp_dir)

        if use_ante_mmpbsa:
            radii = radii_map[str(igb)]
            dr_prmtop = f"{temp_dir}dr.prmtop"
            dl_prmtop = f"{temp_dir}dl.prmtop"
            dc_prmtop = f"{temp_dir}dc.prmtop"
            ante_cmd = f'ante-MMPBSA.py -p {self.prmtop_path} --radii {radii} -s ":WAT,Na+,Cl-" -c {dc_prmtop} -n "{ligand_mask}" -l {dl_prmtop} -r {dr_prmtop}'
            try:
                run(ante_cmd, check=0, text=True, shell=True, capture_output=True)
            except CalledProcessError as err:
                print(err.stdout, err.stderr)
                raise err

        else:
            frag1_path = f"{temp_dir}frag1.pdb"
            frag2_path = f"{temp_dir}frag2.pdb"
            dc_path = f"{temp_dir}dry_complex.pdb"
            # decode ligand mask
            ligand_idx = int(ligand_mask.strip()[1:])
            print(f'working on binding of {ligand_idx}')

            # make new pdb files 
            stru1 = Structure.fromPDB(self.path).find_idx_residue(ligand_idx)
            stru2 = Structure.fromPDB(self.path)
            stru1.build(frag1_path)
            stru2.solvents = []
            stru2.metalatoms = list(filter(lambda x: x.resi_name not in ["Na+", "K+"], stru2.metalatoms))
            stru2.delete_idx_ligand(ligand_idx)
            stru2.build(frag2_path)
            stru3 = Structure.fromPDB(self.path)
            stru3.solvents = []
            stru3.metalatoms = list(filter(lambda x: x.resi_name not in ["Na+", "K+"], stru2.metalatoms))
            stru3.build(dc_path)

            # convert frag pdbs to prmtop files
            dl_prmtop = PDB(frag1_path, wk_dir=self.dir).PDB2FF(
                local_lig = 0,
                igb=igb, ifsolve=0,
                prm_out_path=f'{temp_dir}dl.prmtop', if_prm_only=1)[0]
            dr_prmtop = PDB(frag2_path, wk_dir=self.dir).PDB2FF(
                local_lig = 0,
                igb=igb, ifsolve=0, 
                prm_out_path=f'{temp_dir}dr.prmtop', if_prm_only=1)[0]
            dc_prmtop = PDB(dc_path, wk_dir=self.dir).PDB2FF(
                local_lig = 0,
                igb=igb, ifsolve=0, 
                prm_out_path=f'{temp_dir}dc.prmtop', if_prm_only=1)[0]
            # clean
            if Config.debug < 2:
                run('rm leap.log', check=0, text=True, shell=True, capture_output=True)

        sc_prmtop = type(self).update_radii(
            self.prmtop_path,
            out_path=f'{temp_dir}sc.prmtop', igb=igb)

 
        return dr_prmtop, dl_prmtop, dc_prmtop, sc_prmtop

    def run_mmpbsa(
        self,
        dr_prmtop, dl_prmtop, dc_prmtop, sc_prmtop, traj_file,
        in_file: str = None,
        overwrite_in: str = 0,
        if_cluster_job: bool = True,
        cluster: ClusterInterface = None,
        period: int = 30,
        res_setting: Union[dict, None] = None,
        cluster_debug: bool = 0):
        """mvp function for running mmpbsa"""

        if not in_file:
            in_file = Config.Amber.MMPBSA.build_MMPBSA_in()
        else:
            if os.path.exists(in_file):
                if overwrite_in:
                    Config.Amber.MMPBSA.build_MMPBSA_in(in_file)
            else:
                Config.Amber.MMPBSA.build_MMPBSA_in(in_file)
        temp_dir = f'{self.dir}/temp/'
        mkdir(temp_dir)
        mmpbsa_out_path = f'{temp_dir}mmpbsa.dat'

        if if_cluster_job:
            res_keywords = type(self)._get_default_res_setting_mmpbsa(res_setting)
            cmd = f'{Config.get_PC_cmd(res_keywords["node_cores"])} python2 {Config.Amber.MMPBSA.get_MMPBSA_engine()} -O -i {in_file} -o {mmpbsa_out_path} -sp {sc_prmtop} -cp {dc_prmtop} -rp {dr_prmtop} -lp {dl_prmtop} -y {traj_file}'
            mmpbsa_job = job_manager.ClusterJob.config_job(
                commands = cmd,
                cluster = cluster,
                env_settings = cluster.AMBER_ENV['CPU'],
                res_keywords = res_keywords,
                sub_dir = './', # because path are relative
                sub_script_path = f'{temp_dir}/submit_MMPBSA.cmd')
            mmpbsa_job.submit()
            mmpbsa_job.wait_to_end(period=period)
        else:
            raise Exception("only support cluster job mode right now")
        # clean
        if Config.debug < 2:
            run('rm reference.frc _MMPBSA_info', check=0, text=True, shell=True, capture_output=True)
            os.remove("./tmp/MMPBSA.in")
            os.remove(mmpbsa_job.sub_script_path)

        return mmpbsa_out_path

    @staticmethod
    def _get_default_res_setting_mmpbsa(res_setting):
        """TODO combine those repeating functions"""
        if res_setting is None:
            res_setting = dict()
        if isinstance(res_setting, dict):
            res_setting_holder = copy.deepcopy(Config.Amber.MMPBSA.RES)
            # replace assigned key in default dict
            for k, v in res_setting.items():
                res_setting_holder[k] = v
            return res_setting_holder
        if isinstance(res_setting, str):
            return res_setting

    @staticmethod
    def extract_mmpbsa_out(mmpbsa_out_file: str) -> Dict[str, pd.DataFrame]:
        """mvp function extracting mmpbsa out file"""
        gb_pattern = r"GENERALIZED BORN:(?:.|\n)+?Differences \(Complex - Receptor - Ligand\):\n((?:.|\n)+?)-------------------------------------------------------------------------------\n-------------------------------------------------------------------------------"
        pb_pattern = r"POISSON BOLTZMANN:(?:.|\n)+?Differences \(Complex - Receptor - Ligand\):\n((?:.|\n)+?)-------------------------------------------------------------------------------\n-------------------------------------------------------------------------------"
        with open(mmpbsa_out_file) as f:
            f_str = f.read()
            gb_result_tb = re.search(gb_pattern, f_str).group(1).strip()
            pb_result_tb = re.search(pb_pattern, f_str).group(1).strip()
            gb_result = PDB._extract_pb_gb_table(gb_result_tb)
            pb_result = PDB._extract_pb_gb_table(pb_result_tb)
        return {"pb":pb_result, "gb":gb_result}
    
    @staticmethod
    def _extract_pb_gb_table(table_str: str) -> pd.DataFrame:
        """the table looks like this:
        Energy Component            Average              Std. Dev.   Std. Err. of Mean
        -------------------------------------------------------------------------------
        VDWAALS                    -20.3834                3.1308              0.3131
        EEL                        -19.1765                7.2973              0.7297
        EGB                         20.8964                4.2728              0.4273
        ESURF                       -3.8046                0.1256              0.0126

        DELTA G gas                -39.5599                8.4658              0.8466
        DELTA G solv                17.0918                4.2257              0.4226

        DELTA TOTAL                -22.4680                5.2186              0.5219
        Returns: (like this) pd.DataFrame
                          mean      sd     sem
            VDWAALS      -20.3834  3.1308  0.3131
            EEL          -19.1765  7.2973  0.7297
            EGB           20.8964  4.2728  0.4273
            ESURF         -3.8046  0.1256  0.0126
            DELTA G gas  -39.5599  8.4658  0.8466
            DELTA G solv  17.0918  4.2257  0.4226
            DELTA TOTAL  -22.4680  5.2186  0.5219"""

        data_pattern = r"([A-z]+(?: [A-z]+)*)( *[0-9\-\.]+)( *[0-9\-\.]+)( *[0-9\-\.]+)"
        data_lines = re.findall(data_pattern, table_str)
        data_dict = {}
        for line in data_lines:
            data_dict[line[0].strip()] = [float(line[1].strip()), float(line[2].strip()), float(line[3].strip())]
        
        result_df = pd.DataFrame.from_dict(
            data_dict, 
            orient="index",
            columns=["mean", "sd", "sem"])
        
        return result_df
        
    @staticmethod
    def update_radii(prmtop_path, out_path, igb: int = 5):
        '''
        Update Radii for the MMGBSA calculation
        '''
        # get radii
        radii = radii_map[str(igb)]

        temp_dir = './parmed_tmp/'
        mkdir(temp_dir)

        new_prmtop_path = out_path
        # change Radii
        with open(f'{temp_dir}/parmed.in','w') as of:
            of.write('changeRadii '+radii+line_feed)
            of.write('parmout '+new_prmtop_path+line_feed)
        try:
            run(f'parmed -O -p {prmtop_path} -i {temp_dir}/parmed.in', 
                check=True, text=True, shell=True, capture_output=True)
        except CalledProcessError as err:
            print(err.stdout, err.stderr)
            raise err

        # clean
        if Config.debug < 2:
            rmtree(temp_dir)

        return out_path


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
    - Test result: `cat 1NVG.pdb |grep "^ATOM\\|^HETATM\\|^TER|^END" > 1NVG-grep.pdb` will be fine
    - WARNINGs
    '''
    pass

