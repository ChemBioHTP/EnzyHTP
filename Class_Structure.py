import numpy as np
import os, re
from math import ceil
from Class_line import PDB_line
from Class_Conf import Config
from helper import Child, get_center, get_distance, line_feed, mkdir
from AmberMaps import *
try:
    import openbabel
    import openbabel.pybel as pybel
except ImportError:
    raise ImportError('OpenBabel not installed.')
__doc__='''
This module extract and operate structural infomation from PDB
# will replace some local function in PDB class in the future.
-------------------------------------------------------------------------------------
Class Structure
-------------------------------------------------------------------------------------
Class Chain
-------------------------------------------------------------------------------------
Class Residue
-------------------------------------------------------------------------------------
Class Atom
-------------------------------------------------------------------------------------
Class Metalatom
-------------------------------------------------------------------------------------
Class Ligand
-------------------------------------------------------------------------------------
Class Solvent
===============
'''

class Structure():
    '''
    initilize from
    PDB:        Structure.fromPDB(input_obj, input_type='path' or 'file' or 'file_str')
    raw data:   Structure(chains=None, metalatoms=None, ligands=None)
    ------------
    chains = [chain_obj, ...]
    metalatoms = [metalatom_obj, ...]
    ligands = [ligand, ...]
    solvents = [solvent, ...]
    # maybe add solvent in the future. mimic the metal treatment
    ------------
    METHOD
    ------------
    if_metalatom

    # Move in in the future
    if_art_resi
    if_ligand
    if_complete
    ligand

    get_metal_center
    get_art_resi
    add
    sort
    build

    protonation_metal_fix

    get_all_protein_atom

    ---------
    Special Method
    ---------
    __len__
        len(obj) = len(obj.child_list)
    TODO support better way to operate chains metelatoms ligands (e.g.: search by name, ... instead of just use the list id) (*The list id will have problem since the ligand, metal, and solvent are read in a reversed manner.)
    '''

    '''
    ====
    init
    ====
    '''

    def __init__(self, chains=[], metalatoms=[], ligands=[], solvents=[], name = None):
        '''
        Common part of init methods: direct from data objects
        '''
        if len(chains) == 0 and len(metalatoms) == 0 and len(ligands) == 0 and len(solvents) == 0:
            raise ValueError('need at least one input')
            
        self.chains = []
        self.metalatoms = []
        self.ligands = []
        self.metal_centers = []
        self.solvents = []

        # Add parent pointer and combine into a whole
        for chain in chains:
            chain.set_parent(self)
            self.chains.append(chain)
        for metalatom in metalatoms:
            metalatom.set_parent(self)
            self.metalatoms.append(metalatom)
        for ligand in ligands:
            ligand.set_parent(self)
            self.ligands.append(ligand)
        for solvent in solvents:
            solvent.set_parent(self)
            self.solvents.append(solvent)
        self.name = name

    @classmethod
    def fromPDB(cls, input_obj, input_type='path', input_name = None, ligand_list = None):
        '''
        extract the structure from PDB path. Capable with raw experimental and Amber format
        ---------
        input = path (or file or file_str)
            split the file_str to chain and init each chain
        ligand_list: ['NAME',...]
            User specific ligand names. Only extract these if provided. 
        ---------
        Target:
        - structure(w/name)  - chain - residue - atom
                            |- metalatom(atom)
                            |- ligand(residue)
                            |- solvent(residue)
        - ... (add upon usage)
        ''' 

        # adapt general input // converge to file_str
        if input_type == 'path':
            f = open(input_obj)
            file_str = f.read()
            f.close()
        if input_type == 'file':
            file_str = input_obj.read()
        if input_type == 'file_str':
            file_str = input_obj
        
        raw_chains = []
        # get raw chains
        chains_str = file_str.split(line_feed+'TER') # Note LF is required
        for index, chain_str in enumerate(chains_str):
            if chain_str.strip() != 'END' and chain_str.strip() != '':
                Chain_index = chr(65+index) # Covert to ABC using ACSII mapping
                # Generate chains
                raw_chains.append(Chain.fromPDB(chain_str, Chain_index))
        # clean chains
        # clean metals
        raw_chains_woM, metalatoms = cls._get_metalatoms(raw_chains, method='1')
        # clean ligands
        raw_chains_woM_woL, ligands = cls._get_ligands(raw_chains_woM, ligand_list=ligand_list)
        # clean solvent
        raw_chains_woM_woL_woS, solvents = cls._get_solvents(raw_chains_woM_woL)

        ####### debug ##########
        if Config.debug > 1:
            for chain in raw_chains_woM_woL_woS:
                print('Structure.fromPDB: final chain sequence: '+chain.id, chain.get_chain_seq(Oneletter=1))
            for metal in metalatoms:
                print('Structure.fromPDB: final metal recorded '+metal.name)
            for ligand in ligands:
                print('Structure.fromPDB: final ligand recorded ' + ligand.name)

        return cls(raw_chains_woM_woL_woS, metalatoms, ligands, solvents, input_name)


    @classmethod
    def _get_metalatoms(cls, raw_chains, method='1'):
        '''
        get metal from raw chains and clean chains by deleting the metal part
        -----
        Method 1:   Assume metal can be in any chain
                    Assume all resiude have unique index. 
        (Slow but general)
        
        Method 2: Assume metal can only be in a seperate chain
        (fast but limited)
        '''
        metalatoms = []
        if method == '1':
            for chain in raw_chains:
                for i in range(len(chain)-1,-1,-1):
                    # operate in residue level
                    residue = chain[i]
                    if residue.name in Metal_map.keys():
                        # add a logger in the future
                        print('\033[1;34;40mStructure: found metal in raw: '+chain.id+' '+residue.name+' '+str(residue.id)+' \033[0m')
                        metalatoms.append(residue)
                        del chain[i]
        if method == '2':
            # Not finished yet
            pass

        # Break pseudo residues into atoms and convert to Metalatom object 
        holders = []
        for pseudo_resi in metalatoms:
            for metal in pseudo_resi:
                holders.append(Metalatom.fromAtom(metal))
        metalatoms = holders
                
        # clean empty chains
        for i in range(len(raw_chains)-1,-1,-1):
            if len(raw_chains[i]) == 0:
                del raw_chains[i]

        return raw_chains, metalatoms

    @classmethod
    def _get_ligands(cls, raw_chains, ligand_list = None):
        '''
        get ligand from raw chains and clean chains by deleting the ligand part
        -----
        ligand_list: only record user specified ligand if provided
        Method: Assume metal/ligand/solvent can only be in a seperate chain (or it can not be distinguish from artificial residues.)
                - delete names from rd_non_ligand_list
                + get names from ligand_list only if provided.
        '''
        ligands = []
        for chain in raw_chains:
            #determine if a metal/ligand/solvent chain
            if_HET_chain = 1
            for resi in chain:
                if resi.name in Resi_map2:
                   if_HET_chain = 0
                   break
            if if_HET_chain:
                for i in range(len(chain)-1,-1,-1):
                    # operate in residue level
                    residue = chain[i]

                    # User defined ligand
                    if ligand_list is not None:
                        if residue.name in ligand_list:
                            # add a logger in the future
                            print('\033[1;34;40mStructure: found user assigned ligand in raw: '+chain.id+' '+residue.name+' '+str(residue.id)+' \033[0m')
                            ligands.append(residue)
                            del chain[i]
                    else:
                        if residue.name not in rd_solvent_list:
                            if residue.name not in rd_non_ligand_list:
                                print('\033[1;34;40mStructure: found ligand in raw: '+chain.id+' '+residue.name+' '+str(residue.id)+' \033[0m')
                                ligands.append(residue)
                            del chain[i]


        # Convert pseudo residues to Ligand object 
        holders = []
        for pseudo_resi in ligands:
            holders.append(Ligand.fromResidue(pseudo_resi))
        ligands = holders
                
        # clean empty chains
        for i in range(len(raw_chains)-1,-1,-1):
            if len(raw_chains[i]) == 0:
                del raw_chains[i]

        return raw_chains, ligands

    @classmethod
    def _get_solvents(cls, raw_chains):
        '''
        get solvent from raw chains and clean chains by deleting the solvent part
        -----
        Method: Assume metal/ligand/solvent can anywhere. Base on rd_solvent_list
        '''
        solvents = []
        for chain in raw_chains:
            for i in range(len(chain)-1,-1,-1):
                # operate in residue level
                residue = chain[i]
                if residue.name in rd_solvent_list:
                    if Config.debug > 1:
                        print('\033[1;34;40mStructure: found solvent in raw: '+residue.name+' '+str(residue.id)+' \033[0m')
                    solvents.append(residue)
                    del chain[i]

        # Convert pseudo residues to Ligand object 
        holders = []
        for pseudo_resi in solvents:
            holders.append(Solvent.fromResidue(pseudo_resi))
        solvents = holders
                
        # clean empty chains
        for i in range(len(raw_chains)-1,-1,-1):
            if len(raw_chains[i]) == 0:
                del raw_chains[i]

        return raw_chains, solvents

    '''
    ====
    Methods
    ====
    '''
    def get_metal_center(self):
        '''
        Extract metal centers from metalatoms. Judged by the MetalCenter_map
        save to self.metal_centers
        return self.metal_centers
        '''
        self.metal_centers = []
        for metal in self.metalatoms:
            if metal.resi_name in MetalCenter_map:
                self.metal_centers.append(metal)
        return self.metal_centers


    def get_art_resi(self):
        '''
        find art_resi
        '''
        pass


    def add(self, obj, id=None, sort=0):
        '''
        1. judge obj type (go into the list)
        2. assign parent
        3. id
        if sort:
            clean original id (use a place holder to represent last)
        if not None:
            assign id
        if sort and not None:
            mark as id+i
        4. add to corresponding list
                                                    sort
                  |     |               0             |                     1               |
         assigned |  0  |   keep (for direct output)  | clean (for sort)                    |
                  |  1  |  assign (for direct output) | mark  (for relative order in sort ) |
        '''
        # list
        if type(obj) == list:
            
            obj_ele=obj[0]

            if type(obj_ele) != Chain and type(obj_ele) != Metalatom and type(obj_ele) != Ligand and type(obj_ele) != Solvent:
                raise TypeError('structure.Add() method only take Chain / Metalatom / Ligand / Solvent')

            # add parent and clean id (if sort) assign id (if assigned) leave mark if sort and assigned
            #                         sort
            #          |     |    0     |   1   |
            # assigned |  0  |   keep   | clean |
            #          |  1  |  assign  | mark  |
            for i in obj:               
                i.set_parent(self)
                if sort:
                    if id != None:
                        i.id = str(id)+'i' #str mark
                    else:
                        i.id = id #None 
                else:
                    if id != None:
                        i.id=id
            
            if type(obj_ele) == Chain:
                self.chains.extend(obj)
            if type(obj_ele) == Metalatom:
                self.metalatoms.extend(obj)
            if type(obj_ele) == Ligand:
                self.ligands.extend(obj)
            if type(obj_ele) == Solvent:
                self.solvents.extend(obj)

        # single building block
        else:
            if type(obj) != Chain and type(obj) != Metalatom and type(obj) != Ligand and type(obj) != Solvent:
                raise TypeError('structure.Add() method only take Chain / Metalatom / Ligand / Solvent')
            
            obj.set_parent(self)
            if sort:
                if id != None:
                    obj.id = str(id)+'i' #str mark
                else:
                    obj.id = id #None 
            else:
                if id != None:
                    obj.id=id

            if type(obj) == Chain:
                self.chains.append(obj)
            if type(obj) == Metalatom:
                self.metalatoms.append(obj)   
            if type(obj) == Ligand:
                self.ligands.append(obj)
            if type(obj) == Solvent:
                self.solvents.append(obj)
            
        if sort:
            self.sort()
            

    def sort(self, if_local = 0):
        '''
        assign index according to current items
        chain.id
        resi.id
        atom.id
        -----------
        Chain/Residue level: 
            Base on the order of the old obj.id 
            and potential insert mark from add (higher than same number without the mark)
            *if added object has same id and is not assigned with a insert mark -- place after a original one.
        Atom level:
            base on the parent order (parent.id):
            chains -> metalatoms -> ligands
            residue.id within each above.
            list order within each residues.
        '''
        if if_local:
            # sort chain order
            self.chains.sort(key=lambda chain: chain.id)
            # rename each chain
            for index, chain in enumerate(self.chains):
                chain.id = chr(65+index) # Covert to ABC using ACSII mapping
                # sort each chain
                chain.sort()

            # sort ligand // Do I really need?
            for ligand in self.ligands:
                ligand.sort() #Do nothing
        else:
            r_id = 0
            a_id = 0
            for chain in self.chains:
                for res in chain:
                    r_id +=1
                    res.id = r_id
                    for atom in res:
                        a_id +=1
                        atom.id = a_id
            for metal in self.metalatoms:
                a_id +=1
                metal.id = a_id
            for lig in self.ligands:
                r_id +=1
                lig.id = r_id
                for atom in lig:
                    a_id +=1
                    atom.id = a_id
            for sol in self.solvents:
                r_id += 1
                sol.id = r_id
                for atom in sol:
                    a_id += 1
                    atom.id = a_id


    def build(self, path, ff='AMBER', forcefield='ff14SB', keep_id = 0):
        '''
        build PDB after the change based on the chosen format and forcefield
        - line based on atom and contain chain index and residue index
        ----------------------------
        ff = 
        AMBER (standard amber format: from tleap examples)
            - resi and atom indexes start from 1 and DO NOT reset reaching a new chain. 
            - use atom and residue names from amber force field.
            - place metal, ligand, solvent in seperate chains (seperate with TER)
            - ligand -> metal -> solvent order
            * do not sort atomic order in a residue like tleap does.
        '''
        with open(path, 'w') as of:
            if ff == 'AMBER':
                if not keep_id:
                    a_id = 0
                    r_id = 0
                    for chain in self.chains:
                        #write chain
                        for resi in chain:
                            r_id = r_id+1
                            for atom in resi:
                                a_id = a_id + 1 #current line index
                                line = atom.build(a_id= a_id, r_id = r_id, ff=ff, forcefield=forcefield)
                                of.write(line)
                        #write TER after each chain
                        of.write('TER'+line_feed)

                    c_id = chr(len(self.chains)+64)

                    for ligand in self.ligands:
                        r_id = r_id + 1
                        c_id = chr(ord(c_id)+1)

                        for atom in ligand:
                            a_id = a_id + 1
                            line = atom.build(a_id= a_id, r_id = r_id, c_id = c_id, ff=ff, forcefield=forcefield)
                            of.write(line)
                        of.write('TER'+line_feed)

                    for metal in self.metalatoms:
                        a_id = a_id + 1
                        r_id = r_id + 1
                        c_id = chr(ord(c_id)+1)

                        line = metal.build(a_id= a_id, r_id = r_id, c_id = c_id, ff=ff, forcefield=forcefield)
                        of.write(line)
                        of.write('TER'+line_feed)

                    if len(self.solvents) != 0:
                        c_id = chr(ord(c_id)+1) # same chain_id for all solvent
                        for solvent in self.solvents:
                            r_id = r_id + 1
                            for atom in solvent:
                                a_id = a_id + 1
                                line = atom.build(a_id= a_id, r_id = r_id, c_id = c_id, ff=ff, forcefield=forcefield)
                                of.write(line)
                        of.write('TER'+line_feed)

                else:
                    for chain in self.chains:
                        #write chain
                        for resi in chain:
                            for atom in resi:
                                line = atom.build(ff=ff, forcefield=forcefield)
                                of.write(line)
                        #write TER after each chain
                        of.write('TER'+line_feed)

                    c_id = chr(len(self.chains)+64)

                    for ligand in self.ligands:
                        c_id = chr(ord(c_id)+1)
                        for atom in ligand:
                            line = atom.build(c_id = c_id, ff=ff, forcefield=forcefield)
                            of.write(line)
                        of.write('TER'+line_feed)

                    for metal in self.metalatoms:
                        c_id = chr(ord(c_id)+1)
                        line = metal.build(c_id = c_id, ff=ff, forcefield=forcefield)
                        of.write(line)
                        of.write('TER'+line_feed)

                    if len(self.solvents) != 0:
                        c_id = chr(ord(c_id)+1) # chain_id for all solvent
                        for solvent in self.solvents:
                            for atom in solvent:
                                line = atom.build(c_id = c_id, ff=ff, forcefield=forcefield)
                                of.write(line)
                        of.write('TER'+line_feed)


            if ff == 'XXX':
                #place holder
                pass
                
            of.write('END'+line_feed)


    def build_ligands(self, dir, ft='PDB', ifcharge=0 ,c_method='PYBEL', ph=7.0, ifname=0, ifunique=0):
        '''
        build files for every ligand in self.ligands
        -------
        dir      : output dir. (e.g. File path for ligand i is $dir/ligand_i.pdb)
        ft       : file type / now support: PDB(default)
        ifcharge : if calculate net charge info. (do not before add H)
        c_method : method determining the net charge (default: PYBEL)
        ph       : pH value used for determine the net charge
        ifname   : export residue name if 1 (default: 0)
        ifunique : 1: only build one ligand if there's multiple same ones. 0: build every ligand
        '''
        out_ligs = []

        l_id = 0
        lig_list = self.ligands
        # Only build ligand with same name once
        if ifunique:
            lig_list_uni = {}
            for lig in lig_list:
                lig_list_uni[lig.name] = lig
            lig_list = lig_list_uni.values()

        for lig in lig_list:
            l_id = l_id + 1 # current ligand id
            # make output path
            if dir[-1] == '/':
                dir = dir[:-1]
            if ifunique:
                out_path = dir+'/ligand_'+lig.name+'.pdb'
            else:
                out_path = dir+'/ligand_'+str(l_id)+'_'+lig.name+'.pdb'
            # write
            lig.build(out_path, ft=ft)
            # net charge
            net_charge=None
            if ifcharge:
                if lig.net_charge != None:
                    net_charge = lig.net_charge
                else:
                    net_charge = lig.get_net_charge(method=c_method, ph=ph, o_dir=dir)

            # record
            if ifname:
                out_ligs.append((out_path, net_charge, lig.name))
            else:
                out_ligs.append((out_path, net_charge))

        return out_ligs


    def build_protein(self, dir, ft='PDB'):
        '''
        build only protein and output under the dir
        -------
        dir: out put dir ($dir/protein.pdb)
        ft : file type / now support: PDB(default) 
        '''
        # make path
        if dir[-1] == '/':
            dir = dir[:-1]
        out_path = dir+'/protein.pdb'

        #write
        if ft == 'PDB':
            with open(out_path,'w') as of:
                a_id = 0
                r_id = 0
                for chain in self.chains:
                    #write chain
                    for resi in chain:
                        r_id = r_id+1
                        for atom in resi:
                            a_id = a_id + 1 #current line index
                            line = atom.build(a_id= a_id, r_id = r_id)
                            of.write(line)
                    #write TER after each chain
                    of.write('TER'+line_feed)
                of.write('END'+line_feed)
        else:
            raise Exception('Support only PDB output now.')

        return out_path
    
    
    def build_metalcenters(self, dir, ft='PDB'):
        '''
        build metalcenters only. Use for MCPB parameterization. Deal with donor residue with different protonation states.        ----------
        TODO
        '''
        out_paths = []
        return out_paths


    def get_connect(self, metal_fix = 1, ligand_fix = 1, prepi_path=None):
        '''
        get connectivity
        -----------------
        TREATMENT
        chain: based on connectivity map of each atom in each residue
        metalatom:  fix1: treat as isolated atom
                    fix2: connect to donor atom (MCPB?)
        ligand: fix1: use antechamber generated prepin file to get connectivity.
                      according to https://ambermd.org/doc/prep.html the coordniate line will always start at the 11th line after 3 DUMM.
        '''
        # san check
        if ligand_fix == 1 and prepi_path == None:
            raise Exception('Ligand fix 1 requires prepin_path.')
        # chain part
        for chain in self.chains:
            for res in chain:
                for atom in res:
                    atom.get_connect()
        for sol in self.solvents:
            for atom in sol:
                atom.get_connect()
        # metal
        for metal in self.metalatoms:
            metal.connect = []
        if metal_fix == 1:
            pass
        if metal_fix == 2:
            raise Exception('TODO: Still working on 2 right now')
        
        # ligand
        # init
        for lig in self.ligands:
            for atom in lig:
                atom.connect = []
        # fix 1
        if ligand_fix == 1:
            for lig in self.ligands:
                # read prepin for each ligand
                with open(prepi_path[lig.name]) as f:
                    line_id = 0
                    if_loop = 0
                    for line in f:
                        line_id += 1
                        if line.strip() == '':
                            if if_loop == 1:
                                # switch off loop and break if first blank after LOOP encountered
                                if_loop = 0
                                break
                            continue
                        if if_loop:
                            lp = line.strip().split()
                            lig._find_atom_name(lp[0]).connect.append(lig._find_atom_name(lp[1]))
                            continue
                        # loop connect starts at LOOP
                        if line.strip() == 'LOOP':
                            if_loop = 1
                            continue
                        # coord starts at 11th
                        if line_id >= 11:
                            lp = line.strip().split()
                            atom_id = int(lp[0])-3
                            atom_cnt = int(lp[4])-3
                            if atom_cnt != 0:
                                lig[atom_id-1].connect.append(lig[atom_cnt-1])
                                lig[atom_cnt-1].connect.append(lig[atom_id-1])
                            

    def get_connectivty_table(self, ff='GAUSSIAN', metal_fix = 1, ligand_fix = 1, prepi_path=None):
        '''
        get connectivity table with atom index based on 'ff' settings:
        ff = GAUSSIAN  -- continuous atom index start from 1, do not seperate by chain
        -------------------
        TREATMENT
        Use original atom.id.
            chain: based on connectivity map of each atom in each residue
            metalatom:  fix1: treat as isolated atom
                        fix2: connect to donor atom (MCPB?)
            ligand: fix1: use antechamber generated prepin file to get connectivity.
        Use 1.0 for all connection.
            Amber force field in gaussian do not account in bond order. (Only UFF does.)
            Note that bond order less than 0.1 do not count in MM but only in opt redundant coordinate.
        '''
        connectivty_table = ''
        # get connect for every atom in stru
        self.get_connect(metal_fix, ligand_fix, prepi_path)

        # write str in order
        # Note: Only write the connected atom with larger id
        a_id = 0
        for chain in self.chains:
            for res in chain:
                for atom in res:
                    a_id += 1
                    cnt_line = ' '+str(atom.id)
                    # san check
                    if atom.id != a_id:
                        raise Exception('atom id error.')
                    for cnt_atom in atom.connect:
                        if cnt_atom.id > atom.id:
                            cnt_line += ' '+str(cnt_atom.id)+' '+'1.0'
                    connectivty_table += cnt_line+line_feed
                        
        for lig in self.ligands:
            for atom in lig:
                a_id += 1
                cnt_line = ' '+str(atom.id)
                # san check
                if atom.id != a_id:
                    raise Exception('atom id error.')
                for cnt_atom in atom.connect:
                    if cnt_atom.id > atom.id:
                        cnt_line += ' '+str(cnt_atom.id)+' '+'1.0'
                connectivty_table += cnt_line+line_feed

        for atom in self.metalatoms:
            a_id += 1
            cnt_line = ' '+str(atom.id)
            # san check
            if atom.id != a_id:
                raise Exception('atom id error.')
            for cnt_atom in atom.connect:
                if cnt_atom.id > atom.id:
                    cnt_line += ' '+str(cnt_atom.id)+' '+'1.0'
            connectivty_table += cnt_line+line_feed

        for sol in self.solvents:
            for atom in sol:
                a_id += 1
                cnt_line = ' '+str(atom.id)
                # san check
                if atom.id != a_id:
                    raise Exception('atom id error.')
                for cnt_atom in atom.connect:
                    if cnt_atom.id > atom.id:
                        cnt_line += ' '+str(cnt_atom.id)+' '+'1.0'
                connectivty_table += cnt_line+line_feed

        return connectivty_table


    def protonation_metal_fix(self, Fix):
        '''
        return a bool: if there's any metal center
        '''
        # try once if not exist
        if self.metal_centers == []:
            self.get_metal_center()
        if self.metal_centers == []:
            print('No metal center is found. Exit Fix.')
            return False

        # start fix
        # get donor atoms and residues
        for metal in self.metal_centers:
            metal.get_donor_residue(method = 'INC')

            if Fix == 1:
                metal._metal_fix_1()
            
            if Fix == 2:
                metal._metal_fix_2()

            if Fix == 3:
                metal._metal_fix_3()
        return True


    def get_all_protein_atom(self):
        '''
        get a list of all protein atoms
        return all_P_atoms 
        '''
        all_P_atoms = []
        for chain in self.chains:
            for residue in chain:
                all_P_atoms.extend(residue.atoms)
        return all_P_atoms


    def get_all_residue_unit(self, ifsolvent=0):
        all_r_list = []
        for chain in self.chains:
            for resi in chain:
                all_r_list.append(resi)

        for resi in self.ligands:
            all_r_list.append(resi)

        for resi in self.metalatoms:
            all_r_list.append(resi)

        if ifsolvent:
            for resi in self.solvents:
                all_r_list.append(resi)

        return all_r_list


    def get_residue(self, id):
        '''
        re-search residue id with all residues count togethor from 1.
        ----------
        return a residue object
        '''
        all_resi = self.get_all_residue_unit()
        for resi in all_resi:
            if resi.id == int(id):
                return resi


    def get_atom_id(self):
        '''
        return a list of id of all atoms in the structure
        '''
        atom_id_list = []
        for chain in self.chains:
            for res in chain:
                for atom in res:
                    if atom.id == None:
                        raise Exception('Detected None in chain '+str(chain.id)+str(res.id)+' '+atom.name)
                    atom_id_list.append(atom.id)
        
        for metal in self.metalatoms:
            if metal.id == None:
                raise Exception('Detected None in metal '+ metal.name)
            atom_id_list.append(metal.id)
        
        for lig in self.ligands:
            for atom in lig:
                if atom.id == None:
                    raise Exception('Detected None in ligands', res.id, atom.name)
                atom_id_list.append(atom.id)

        for sol in self.solvents:
            for atom in sol:
                if atom.id == None:
                    raise Exception('Detected None in solvent', res.id, atom.name)
                atom_id_list.append(atom.id)

        return atom_id_list


    def get_atom_charge(self, prmtop_path):
        '''
        requires generate the stru using !SAME! PDB as one that generate the prmtop. 
        '''
        pass   


    def get_atom_type(self, prmtop_path):
        '''
        requires generate the stru using !SAME! PDB as one that generate the prmtop. 
        '''
        # get type list
        with open(prmtop_path) as f:
            type_list=[]
            line_index=0

            for line in f:
                
                line_index=line_index+1 #current line
                
                if line.strip() == r'%FLAG POINTERS':
                    format_flag=line_index
                if line.strip() == r'%FLAG AMBER_ATOM_TYPE':
                    type_flag=line_index

                if 'format_flag' in dir():
                    if line_index == format_flag+2:
                        N_atom=int(line.split()[0])
                        del format_flag
                        
                if 'type_flag' in dir():
                    if line_index >= type_flag+2 and line_index <= type_flag+1+ceil(N_atom/5):
                        for i in line.strip().split():
                            type_list.append(i)
        # assign type to atom
        for chain in self.chains:
            for res in chain:
                for atom in res:
                    atom.type = type_list[atom.id-1]
        for atom in self.metalatoms:
            if atom.id == None:
                raise Exception('Detected None in metal '+ atom.name)
            atom.type = type_list[atom.id-1]
        for lig in self.ligands:
            for atom in lig:
                if atom.id == None:
                    raise Exception('Detected None in ligands', res.id, atom.name)
                atom.type = type_list[atom.id-1]
        for sol in self.solvents:
            for atom in sol:
                if atom.id == None:
                    raise Exception('Detected None in solvent', res.id, atom.name)
                atom.type = type_list[atom.id-1]
        
        
    def get_sele_list(self, atom_mask, fix_end='H', prepi_path=None):
        '''
        interface with class ONIOM_Frame. Generate a list for sele build. Make sure use same pdb as the one generate the frame.
        ------------
        resi_list: selected residue list
        atom_mask: atom selection with the standard grammer of Amber (incomplete)
        fix_end: fix valence of the cut bond. (default: H)
                - H: add H to where the original connecting atom is.
                    special fix for classical case:
                    "sele by residue" (cut N-C) -- adjust dihedral for added H on N.
                = Interface with write_sele_lines:
                    add {fix_flag+element_mark: coord} in sele_lines 
        ------------
        return a sele list:
        - Fixing atoms are labeled as qm_atom_id-qm_atom_cnt_id-distance
        - backbone atoms are marked as b at the end (for qmcluster charge calculation)
        - other atoms use _ as place holder
        return a sele map:
        (PDB atom id -> QM atom id)
        '''
        sele_lines = {}
        #decode atom_mask (maybe in helper later) TODO
        resi_list = atom_mask[1:].strip().split(',')
        all_resi_list = self.get_all_residue_unit()

        # decode and get obj
        sele_stru_objs=[]
        for resi in resi_list:
            chain_id = re.match('[A-Z]',resi)
            resi_id = int(re.match('[0-9]+',resi).group(0))
            if chain_id == None:
                for resi in all_resi_list:
                    if resi_id == resi.id:
                        resi_obj = resi
            else:
                chain_id = chain_id.group(0)
                resi_obj = self.chains[int(chain_id)-65]._find_resi_id(resi_id)

            sele_stru_objs.append(resi_obj)

        # combine the sele
        sele_atoms = []
        for obj in sele_stru_objs:
            for atom in obj:
                sele_atoms.append(atom)

        if fix_end != None:
            self.get_connect(prepi_path=prepi_path)

        # operate on the sele objs
        for atom in sele_atoms:
            # add current atom
            atom.get_ele()
            if type(atom.parent) != Ligand:
                if atom.name in ['C','CA','O','N','H','HA']:
                    sele_lines[str(atom.id)+'b'] = atom.ele
                else:
                    sele_lines[str(atom.id)+'_'] = atom.ele
            else:
                sele_lines[str(atom.id)+'_'] = atom.ele
            
            if fix_end != None:
                # search for cut bond
                for cnt_atom in atom.connect:
                    if not cnt_atom in sele_atoms:
                        if fix_end == 'H':
                            d_XH = X_H_bond_length[atom.name]
                            label = '-'.join((str(atom.id),str(cnt_atom.id),str(d_XH)))
                            fix_atom = 'H'
                        if fix_end == 'Me':
                            #TODO
                            pass
                        # write to sele_lines
                        sele_lines[label] = fix_atom
        # make sele_map (PDB atom id -> QM atom id)
        sele_map = {}
        for i, key in enumerate(sele_lines.keys()):
            if key[-1] not in '1234567890':
                key = key[:-1]
            sele_map[key] = i+1

        if Config.debug >= 1:
            print('Selected QM cluster atoms: ')
            print(sele_lines)

        return sele_lines, sele_map


    def get_resi_dist(self, r1, r2, method='mass_center'):
        '''
        r1: residue 1. (residue_obj)
        r2: residue 2. (residue_obj)
        method: The method to narrow the residue down to a point.
                - mass_center
                - ...
        '''
        if method == 'mass_center':
            p1 = r1.get_mass_center()
            p2 = r2.get_mass_center()

        D = get_distance(p1, p2)

        return D



    '''
    ====
    Special Method
    ====
    '''

    def __len__(self):
        '''
        len(obj) = len(obj.child_list)
        '''
        return len(self.chains)+len(self.metalatoms)+len(self.ligands)+len(self.solvents)


    def __getitem__(self, i):
        '''
        pop four elements with fixed order
        -----------------
        capable with future update
        '''

        if i == 0:
            return self.chains
        if i == 1:
            return self.ligands
        if i == 2:
            return self.metalatoms
        if i == 3:
            return self.solvents
        if i > 3:
            raise StopIteration



class Chain(Child):
    '''
    -------------
    initilize from
    PDB:        Chain.fromPDB(chain_input, chain_id, input_type='file_str' or 'file' or 'path')
    raw data:   Chain(residues, chain_id)
    -------------
    id
    parent # the whole structure
    residues = [resi_obj, ...]
    chain_seq = ['resi_name', ..., 'NAN', ..., 'resi_name']
    -------------
    method
    -------------
    set_parent
    get_chain_seq(self)
    _find_resi_name
    -------------
    Special method
    -------------
    __getitem__
        Chain_obj[int]: Chain_obj.residues[int] // (start from 0)
    __getattr__
        Chain_obj.123 = Chain_obj.residues[123-1] // index mimic (start from 1)
        Chain_obj.HIS = Chain_obj.find_resi_name('HIS') // search mimic
    __delitem__
        del obj[int] --> obj.child[int].remove() // delete by index (start from 0)
        del obj[str] --> obj._del_child_name() // delete by value
        del obj[child] --> obj.child_list.remove(child) // delete by value
    __len__
        len(obj) = len(obj.child_list)
    '''

    '''
    ====
    init
    ====
    '''
    def __init__(self, residues, chain_id: str, parent=None):
        '''
        Common part of init methods: direct from data objects
        No parent by default. Add parent by action
        '''
        #set parent to None
        Child.__init__(self)
        # add parent if provided
        if parent != None:
            self.set_parent(parent)
        #adapt some children
        self.residues = []
        for i in residues:
            i.set_parent(self)
            self.residues.append(i)
        #set id
        self.id = chain_id

        # init
        self.ifsorted = 0
    
    @classmethod
    def fromPDB(cls, chain_input, chain_id, input_type='file_str'):
        '''
        generate chain from PDB. Capable with raw experimental and Amber format. Only read 'ATOM' and 'HETATM' lines.
        ---------
        chain_input = file_str (or path or file)
        chain_id : str
        split the file_str to residues and init each residue
        ''' 

        # adapt general input // converge to file_str
        if input_type == 'path':
            f = open(chain_input)
            chain_str = f.read().strip()
            f.close()
        if input_type == 'file':
            chain_str = chain_input.read().strip()
        if input_type == 'file_str':
            chain_str = chain_input.strip()

        # chain residues
        residues = []
        resi_lines = [] # data holder
        lines = PDB_line.fromlines(chain_str) # Note LF is required
        for pdb_l in lines:
            if pdb_l.line_type == 'ATOM' or pdb_l.line_type == 'HETATM':
                # Deal with the first residue
                if len(residues) == 0 and len(resi_lines) == 0:
                    resi_lines.append(pdb_l)
                    last_resi_index = pdb_l.resi_id
                    continue

                # find the beginning of a new residue
                if pdb_l.resi_id != last_resi_index:
                    # Store last resi
                    last_resi = Residue.fromPDB(resi_lines, last_resi_index)
                    residues.append(last_resi)
                    # empty the holder for current resi
                    resi_lines = []

                resi_lines.append(pdb_l)
                                
                # Update for next loop                
                last_resi_index = pdb_l.resi_id

        # Deal with the last residue (Now add the last one after the loop.)
            # if i == len(lines)-1:
        if resi_lines != []:    # deal with blank chain
            last_resi = Residue.fromPDB(resi_lines, resi_lines[-1].resi_id)
            residues.append(last_resi)
        else:
            if Config.debug >= 1:
                print('Chain.fromPDB: find a empty chain: ' + chain_id)

        return cls(residues, chain_id)

    '''
    ====
    Method
    ====
    '''
    def add(self, obj, id=None, sort=0):
        '''
        1. judge obj type
        2. clean original id
        3. add to corresponding list
        '''
        # list
        if type(obj) == list:
            
            obj_ele=obj[0]

            if type(obj_ele) != Residue:
                raise TypeError('chain.Add() method only take Residue')

            # add parent and clean id (if sort) assign id (if assigned) leave mark if sort and assigned
            for i in obj:               
                i.set_parent(self)
                if sort:
                    if id != None:
                        i.id = str(id)+'i' #str mark
                    else:
                        i.id = id #None 
                else:
                    if id != None:
                        i.id=id
            self.residues.extend(obj)
            

        # single building block
        else:
            if type(obj) != Residue:
                raise TypeError('Chain.Add() method only take Residue')
            
            obj.set_parent(self)
            if sort:
                if id != None:
                    obj.id = str(id)+'i' #str mark
                else:
                    obj.id = id #None 
            else:
                if id != None:
                    obj.id=id
            self.residues.append(obj)

        if sort:
            self.sort()

        

    def sort(self, sort_resi = 1):
        '''
        sort_resi: 1/0 -- if or not sort atoms in each residue. (start from 1 for each residue)

        turn residue index into str and sort with list.sort() -- cope with the insert mark
        * start form 1 in each chain
        * if added object has same id and is not assigned with a insert mark -- place after a original one.
        ----
        assign index according to current items
        resi.id
        atom.id
        '''
        # sort residue order
        for i in self.residues:
            if type(i.id) == str:
                raise Exception('Does not support added residue now: update in the future')
            # TODO
        
        self.residues.sort(key=lambda i: i.id)
        # re-id each residue
        for index, resi in enumerate(self.residues):
            resi.id = index+1
            # sort each residue
            if sort_resi:
                resi.sort()
        
        self.ifsorted = 1
      


    def get_chain_seq(self, Oneletter = 0):
        '''
        get chain sequence from the "residues" list. Store in self.seq
        ----------------------------
        Oneletter
        = 0 // use 3-letter format to represet each residue
        = 1 // use 1-letter format to represet each residue
        ----------------------------
        + Use "NAN"/"-" as a filler to store missing residues (detect internal missing)
        - (WARNING) Require the ligand in seperate chains
        - (WARNING) Only support normal chain! (not ligand or solvent)
        ---------------------
        * Note that the missing sequence infomation will be missing after sort(). So get seq before sort to obtain such info.
        * Need to be update after mutation
        '''
        if self.ifsorted:
            if Config.debug >= 1:
                print('Chain.get_chain_seq: WARNING: missing sequence infomation is missing after sort().')

        chain_seq = []
        for resi in self.residues:
            # first resi
            if len(chain_seq) == 0:
                last_id = resi.id
                chain_seq.append(resi.name)
                continue
            # missing sequence
            missing_length = resi.id - last_id - 1
            if missing_length > 0:
                chain_seq = chain_seq + ['NAN',] * missing_length

            chain_seq.append(resi.name)
            last_id = resi.id

        self.seq = chain_seq
        if Oneletter == 1:
            self._get_Oneletter()
            return self.seq_one
        else:
            return self.seq


    def if_art_resi(self):
        '''
        TODO
        '''
        pass


    def _get_Oneletter(self):
        '''
        (Used internally) convert sequences in self.sequence to oneletter-based str
        - The 'NAN' is convert to '-'
        - Present unnature residue as full 3-letter name
        save to self.seq_one
        '''
        if len(self.seq) == 0:
            raise IndexError("The self.sequence should be obtained first")
        
        seq_one = ''
        for name in self.seq:
            if name == 'NAN':
                seq_one=seq_one+'-'
            else:
                if name in Resi_map2:
                    seq_one = seq_one + Resi_map2[name]
                else:
                    seq_one = seq_one + ' '+name+' '

        self.seq_one = seq_one
                
            

    def _find_resi_name(self, name: str):
        '''
        find residues according to the name
        return a list of found residues
        ''' 
        out_list = []
        for resi in self.residues:
            if resi.name == name:
                out_list.append(resi)
        return out_list

    def _find_resi_id(self, id: int):
        '''
        find residues according to the id
        return a found residue
        ''' 
        out_list = []
        for resi in self.residues:
            if resi.id == id:
                out_list.append(resi)
        if len(out_list) > 1:
            print('\033[32;0mShould there be same residue id in chain +'+self.name+str(self.id)+'?+\033[0m')
            raise Exception
        return out_list[0]
    
    def _del_resi_name(self, name: str):
        '''
        find residues according to the name
        delete found residues
        ''' 
        for i in range(len(self.residues)-1,-1,-1):
            if self.residues[i].name == name:
                del self.residues[i]


    '''
    ====
    Special Method
    ====
    '''
    
    def __getitem__(self, key: int):
        '''
        Chain_obj[int]: Chain_obj.residues[int]
        -----
        use residue index within the chain, start from 0
        '''
        return self.residues[key]
    
    def __getattr__(self, key):
        '''
        Use i to indicate id (! by this, do not support searching any residue name start with i)
        Chain_obj.i123 = Chain_obj.residues[123-1] // index mimic (start from 1)
        Chain_obj.HIS = Chain_obj.find_resi_name('HIS') // search mimic
        '''
        if key == 'stru':
            return self.parent
        if key[0] == 'i':
            # digit
            key = int(key[1:])
            return self.residues[key-1]
        else:
            # text
            return self._find_resi_name(key)


    def __delitem__(self, key):
        '''
        del obj[int] --> obj.child[int].remove() // delete by index (start from 0)
        del obj[str] --> obj._del_child_name() // delete by value
        del obj[child] --> obj.child_list.remove(child) // delete by value
        '''
        if type(key) == int:
            del self.residues[key]
        if type(key) == str:
            self._del_resi_name(key)
        if type(key) == Residue:
            self.residues.remove(key)

    def __len__(self):
        '''
        len(obj) = len(obj.child_list)
        '''
        return len(self.residues)
        


class Residue(Child):
    '''
    -------------
    initilize from
    PDB:        Residue.fromPDB(resi_input, resi_id=pdb_l.resi_id, input_type='PDB_line' or 'line_str' or 'file' or 'path')
    raw data:   Residue(atoms, resi_id, resi_name)
    -------------
    id
    name
    parent # the whole chain

    atoms = [atom_obj, ...]

    # --after metalatom.get_donor_resi--
    d_atom (donor atom when work as a ligand)
    a_metal (acceptor metal)
    
    #TODO
    if_art_resi
    -------------
    Method
    -------------
    set_parent
    if_art_resi
    deprotonate
    _find_atom_name
    -------------
    __getitem__
        Residue_obj[int]: Residue_obj.residues[int]    
    __getattr__
        Residue_obj.123 = Residue_obj.atoms[123-1] // index mimic (start from 1)
        Residue_obj.CA = Residue_obj.find_atom_name('CA') // search mimic
    __delitem__
        del obj[int] --> obj.child[int].remove() // delete by index (start from 0)
        del obj[str] --> obj._del_child_name(str) // delete by value
        del obj[child] --> obj.child_list.remove(child) // delete by value
    __len__
        len(obj) = len(obj.child_list)
    '''

    '''
    ====
    init
    ====
    '''
    def __init__(self, atoms, resi_id, resi_name, parent=None):
        '''
        Common part of init methods: direct from data objects
        No parent by default. Add parent by action
        '''
        #set parent to None
        Child.__init__(self) 
        # add parent if provided
        if parent != None:
            self.set_parent(parent)
        #adapt some children
        self.atoms = []
        for i in atoms:
            i.set_parent(self)
            self.atoms.append(i)
        #set id
        self.id = resi_id
        self.name = resi_name

        #clean
        self.d_atom = None
    
    @classmethod
    def fromPDB(cls, resi_input, resi_id=None, input_type='PDB_line'):
        '''
        generate resi from PDB. Require 'ATOM' and 'HETATM' lines.
        ---------
        resi_input = PDB_line (or line_str or file or path)
        resi_id : int (use the number in the line by default // support customize)
        Use PDB_line in the list to init each atom
        '''

        # adapt general input // converge to a list of PDB_line (resi_lines)
        if input_type == 'path':
            f = open(resi_input)
            resi_lines = PDB_line.fromlines(f.read())
            f.close()
        if input_type == 'file':
            resi_lines = PDB_line.fromlines(resi_input.read())
        if input_type == 'line_str':
            resi_lines = PDB_line.fromlines(resi_input)
        if input_type == 'PDB_line':
            resi_lines = resi_input
        
        #clean lines
        for i in range(len(resi_lines)-1,-1,-1):
            if resi_lines[i].line_type != 'ATOM' and resi_lines[i].line_type != 'HETATM':
                if Config.debug > 1:
                    print('Residue.fromPDB: delete error line in input.')
                del resi_lines[i]
    
        # Default resi_id
        if resi_id is None:
            resi_id = resi_lines[0].resi_id
        # get name from first line
        resi_name = resi_lines[0].resi_name
        # get child atoms
        atoms = []
        for pdb_l in resi_lines:
            atoms.append(Atom.fromPDB(pdb_l))
        
        return cls(atoms, resi_id, resi_name)
    

    '''
    ====
    Method
    ====
    '''
    def if_art_resi(self):
        pass

    
    def deprotonate(self, T_atom = None, HIP = 'HIE'):
        '''
        check current protonation state.
        deprotonate if applicable (HIP -> HIE by default)
        ---------
        base on T_atom if provided. (refine in the future)
        '''
        if T_atom == None:
            if self.name != 'HIP':
                depro_info = DeProton_map[self.name]
            else:
                if HIP == 'HIE':
                    depro_info = DeProton_map[self.name][0]
                if HIP == 'HID':
                    depro_info = DeProton_map[self.name][1]

            depro_resi = depro_info[0]
            depro_atom = depro_info[1]

            #delete the proton
            del self[depro_atom]

            #change the name
            self.name = depro_resi

        else:
            # assign target atom
            # only affect operations on residues with differences between donor atom and potential deprotonation atom (HIP HIE HID and ARG)
            if self.name == 'HIP':
                if T_atom.name == 'ND1':
                    del self['HD1']
                    self.name = 'HIE'
                    return
                if T_atom.name == 'NE2':
                    del self['HE2']
                    self.name = 'HID'
                    return
                    
            if self.name == 'HIE':
                if T_atom.name == 'ND1':
                    return
                if T_atom.name == 'NE2':
                    del self['HE2']
                    self.name = 'HID'
                    # self.add_H('ND1')
                    #let leap auto complete by now
                    return

            if self.name == 'HID':
                if T_atom.name == 'ND1':
                    del self['HD1']
                    self.name = 'HIE'
                    # self.add_H('NE2')
                    #let leap auto complete by now
                    return
                if T_atom.name == 'NE2':
                    return

            # find the right H to delete for ARG and LYS
            if self.name == 'ARG':
                raise Exception('ARG detected as donor: '+self.chain.id+str(self.id))

            if self.name == 'LYS':
                # the deleted H has to be HZ1 in name
                D1 = np.linalg.norm(np.array(self.HZ1.coord) - np.array(self.a_metal.coord))
                D2 = np.linalg.norm(np.array(self.HZ2.coord) - np.array(self.a_metal.coord))
                D3 = np.linalg.norm(np.array(self.HZ3.coord) - np.array(self.a_metal.coord))
                x = min(D1,D2,D3)
                if x == D1:
                    del self['HZ1']
                    self.name = 'LYN'
                    return
                if x == D2:
                    del self['HZ2']
                    self.HZ1.name = 'HZ2'
                    self.name = 'LYN'
                    return
                if x == D3:
                    del self['HZ3']
                    self.HZ1.name = 'HZ3'
                    self.name = 'LYN'
                    return
                raise Exception

            depro_info = DeProton_map[self.name]
            depro_resi = depro_info[0]
            depro_atom = depro_info[1]
            #delete the proton
            del self[depro_atom]
            #change the name
            self.name = depro_resi

            
    def add_H(self, T_atom):
        '''
        add H on the corresponding T_atom.
        1. make the H
            find H name
            find H coordinate
        2. add H to the residue
            use self.add()
        '''
        pass


    def rot_proton(self, T_atom):
        '''
        rotate the dihedral relate to the target H-atom bond
        -----
        TODO  
        '''
        
        if self.name == 'TRP':
            raise Exception('Error: TRP detected as donor!!!')

        protons = T_atom.get_protons()
        lp_infos = T_atom.get_lp_infos()
        # rotate to lp direction if proton on T_atom
        if len(protons) != 0:
            for proton in protons:
                bond_end1 =  T_atom.get_bond_end_atom()
                bond_end2 =  bond_end1.get_bond_end_atom()
                proton.set_byDihedral(T_atom, bond_end1, bond_end2, value = lp_infos[0]['D'])


    def ifDeProton(self):
        '''
        check if this residue add or minus proton in a pH range of 1-14. (ambiguous protonation state, potential deprotonation)
        -------
        base on self.name. return a bool
        '''
        return self.name in DeProton_map.keys()
    

    def _find_atom_name(self, name: str):
        '''
        find atom according to the name (should find only one atom)
        return the atom (! assume the uniqueness of name)
        ''' 
        out_list = []
        for atom in self.atoms:
            if atom.name == name:
                out_list.append(atom)
        if len(out_list) > 1:
            print('\033[32;0mShould there be same atom name in residue +'+self.name+str(self.id)+'?+\033[0m')
            raise Exception
        else:
            return out_list[0]

    def _del_atom_name(self, name: str):
        '''
        find atoms according to the name
        delete found atoms
        ''' 
        for i in range(len(self.atoms)-1,-1,-1):
            if self.atoms[i].name == name:
                del self.atoms[i]

    def add(self, obj):
        '''
        1. judge obj type
        2. set parent, add to corresponding list
        --------
        Do not take any id: add to the end of the residue atom list by default
        Do not sort: the order of atom is pointless by far.
        '''
        # list
        if type(obj) == list:
            
            obj_ele=obj[0]

            if type(obj_ele) != Atom:
                raise TypeError('residue.Add() method only take Atom')

            obj.set_parent(self)
            self.atoms.extend(obj)
            

        # single building block
        else:
            if type(obj) != Atom:
                raise TypeError('residue.Add() method only take Atom')
            
            obj.set_parent(self)
            self.atoms.append(obj)

    def sort(self):
        '''
        assign atom id based on current order in the list
        * start from 1 in a residue
        '''
        for index, atom in enumerate(self.atoms):
            atom.id = index+1

    def get_mass_center(self):
        '''
        get mass center of current residue
        '''
        M_t = 0
        M_cw = np.array((0.0,0.0,0.0))
        for atom in self:
            atom.get_ele()
            a_mass = Ele_mass_map[atom.ele]
            M_t += a_mass
            M_cw += a_mass * np.array(atom.coord)
        
        M_center = tuple(M_cw/M_t)

        return M_center

    '''
    ====
    Special Method
    ====
    '''

    def __getitem__(self, key: int):
        '''
        Residue_obj[int]: Residue_obj.residues[int]
        -----
        use residue index within the chain, start from 0
        '''
        return self.atoms[key]
    
    def __getattr__(self, key):
        '''
        Residue_obj.i123 = Residue_obj.atoms[123-1] // index mimic (start from 1)
        Residue_obj.CA = Residue_obj.find_atom_name('CA') // search mimic
        '''
        if key == 'chain':
            return self.parent
        # judge if a digit str, since a str will always be passed
        if key[0] == 'i':
            key = int(key[1:])
            return self.atoms[key-1]
        else:
            return self._find_atom_name(key)

    
    def __delitem__(self, key):
        '''
        del obj[int] --> obj.child[int].remove() // delete by index (start from 0)
        del obj[str] --> obj._del_child_name(str) // delete by value
        del obj[child] --> obj.child_list.remove(child) // delete by value
        '''
        if type(key) == int:
            del self.atoms[key]
        if type(key) == str:
            return self._del_atom_name(key)
        if type(key) == Atom:
            self.atoms.remove(key)

    def __len__(self):
        '''
        len(obj) = len(obj.child_list)
        '''
        return len(self.atoms)



class Atom(Child):
    '''
    -------------
    initilize from
    PDB:        Atom.fromPDB(atom_input, input_type='PDB_line' or 'line_str' or 'file' or 'path')
    raw data:   Atom(atom_name, coord)
    -------------
    id
    name
    coord = [x, y, z]
    ele (obtain by method)

    parent # resi
    -------------
    Method
    -------------
    set_parent
    Add_id # use only after the whole structure is constructed

    get_connect
    get_protons
    get_lp_infos
    get_bond_end_atom
    set_byDihedral
    set_byAngle
    set_byBond

    gen_line

    get_around
    get_ele
    -------------
    '''

    def __init__(self, atom_name: str, coord: list, ff: str, atom_id = None, parent = None):
        '''
        Common part of init methods: direct from data objects
        No parent by default. Add parent by action
        '''
        #set parent to None
        Child.__init__(self)
        # add parent if provided
        if parent != None:
            self.set_parent(parent)
        # get data
        self.name = atom_name
        self.coord = coord
        self.id = atom_id
        self.ff = ff

    @classmethod
    def fromPDB(cls, atom_input, input_type='PDB_line', ff = 'Amber'):
        '''
        generate atom from PDB. Require 'ATOM' and 'HETATM' lines.
        ---------
        atom_input = PDB_line (or line_str or file or path)
        '''
        # adapt general input // converge to a PDB_line (atom_line)
        if input_type == 'path':
            f = open(atom_input)
            atom_line = PDB_line(f.read())
            f.close()
        if input_type == 'file':
            atom_line = PDB_line(atom_input.read())
        if input_type == 'line_str':
            atom_line = PDB_line(atom_input)
        if input_type == 'PDB_line':
            atom_line = atom_input
        
        # get data
        atom_name = atom_line.atom_name
        coord = [atom_line.atom_x,atom_line.atom_y,atom_line.atom_z]
        atom_id = atom_line.atom_id

        return cls(atom_name, coord, ff, atom_id=atom_id)


    '''
    ====
    Method
    ====
    '''

    def get_connect(self):
        '''
        find connect atom base on:
        1. AmberMaps.resi_cnt_map
        2. parent residue name
        ------------
        * require standard Amber format (atom name and C/N terminal name)
        save found list of Atom object to self.connect 
        '''
        self.connect = []

        if self.resi.name in rd_solvent_list:
            name_list = resi_cnt_map[self.resi.name][self.name]
        else:
            r = self.resi
            r1 = self.resi.chain.residues[0]
            rm1 = self.resi.chain.residues[-1]

            if r == r1:
                # N terminal
                name_list = resi_nt_cnt_map[self.resi.name][self.name]
            else:
                if r == rm1:
                    # C terminal
                    name_list = resi_ct_cnt_map[self.resi.name][self.name]
                else:
                    name_list = resi_cnt_map[self.resi.name][self.name]

        for name in name_list:
            try:
                if name not in ['-1C', '+1N']:
                    cnt_atom = self.resi._find_atom_name(name)
                if name == '-1C':
                    cnt_resi = self.resi.chain._find_resi_id(self.resi.id-1)
                    cnt_atom = cnt_resi._find_atom_name('C')
                if name == '+1N':
                    cnt_resi = self.resi.chain._find_resi_id(self.resi.id+1)
                    cnt_atom = cnt_resi._find_atom_name('N')                
                self.connect.append(cnt_atom)
            except IndexError:
                if Config.debug >= 1:
                    print('WARNING: '+self.resi.name+str(self.resi.id)+' should have atom: '+name)
                else:
                    pass

    def get_type(self):
        if self.resi.name in rd_solvent_list:
            self.type = G16_label_map[self.parent.name][self.name]
        else:
            r = self.resi
            r1 = self.resi.chain.residues[0]
            rm1 = self.resi.chain.residues[-1]
            if r == r1 and self.name in ['H1','H2','H3']:
                self.type = 'H'
            else:
                if r == rm1 and self.name == 'OXT':
                    self.type = 'O2'
                else:
                    self.type = G16_label_map[self.parent.name][self.name]
        return self.type
        
    def get_pseudo_H_type(self, old_atom, method=1):
        '''
        determine the connecting pseudo H atom type base on the environment for ONIOM input
        A low layer atom will be replaced to a pesudo H for "model-low" layer modeling
        self    : the atom in the "model" layer
        old_atom: the original connecting atom replacing by H. 
        method  : 
        1 ----- : Only consider some main cases:
                    - truncating CA-CB
                    - truncating C-N
                    (abort for other cases)
        -----------
        TARGET
        -----------
        1. provide reasonable parameters according to the atom type
        2. avoid missing parameters (search ff96)
        3. append missing parameters
        '''
        p_H_type = ''
        if (self.name, old_atom.name) not in [('CA','CB'),('C','N'),('CB','CA'),('N','C')]:
            raise Exception('method 1 only support truncations of CA-CB and C-N')
        else:
            if self.name == 'CA':
                #p_H_type = XXX
                pass


        return p_H_type


    def get_protons(self):
        '''
        get connected proton based on the connectivity map
        ---------
        check if connectivity is obtained. get if not.
        '''
        return []

    def get_lp_infos(self):
        '''
        get lone pair based on the connectivity map
        ---------
        check if connectivity is obtained. get if not.
        '''
        return []

    def get_bond_end_atom(self):
        '''
        get the connected atom that closer to CA in the *network*
        -------
        check if connectivity is obtained. get if not.        
        '''
        return None

    def set_byDihedral(self, A2, A3, A4, value):
        '''
        change the atom (self) coordinate by dihedral.
        '''
        pass

    def set_byAngle(self, A2, A3, value):
        '''
        change the atom (self) coordinate by angle.
        '''
        pass

    def set_byBond(self, A2, value):
        '''
        change the atom (self) coordinate by distance.
        '''
        pass


    def build(self, a_id=None, r_id=None, c_id=None,  ff='AMBER', forcefield = 'ff14SB'):
        '''
        generate an output line. End with LF
        return a line str
        use self.id and parent id if not assigned
        -------
        must use after sort!
        '''
        #default
        if a_id == None:
            a_id = self.id
        if r_id == None:
            r_id = self.resi.id
        if c_id == None:
            c_id = self.resi.chain.id
        
        if ff == 'AMBER':
            #build an amber style line here
            l_type = '{:<6}'.format('ATOM')
            a_index = '{:>5d}'.format(a_id)

            if forcefield == 'ff14SB':
                # ff14SB by default
                #atom name TODO when deal with 2-letter element in artifical residue and ligand
                if len(self.name) > 3:
                    a_name = '{:<4}'.format(self.name)
                else:
                    a_name = '{:<3}'.format(self.name)
                    a_name = ' '+a_name
                r_name = '{:>3}'.format(self.resi.name)
            else:
                raise Exception('Only support ff14SB atom/resiude name now')

            c_index = c_id
            r_index = '{:>4d}'.format(r_id)
            x = '{:>8.3f}'.format(self.coord[0])
            y = '{:>8.3f}'.format(self.coord[1])
            z = '{:>8.3f}'.format(self.coord[2])
        
        #example: ATOM   5350  HB2 PRO   347      32.611  15.301  24.034  1.00  0.00
        line = l_type + a_index +' '+a_name + ' ' + r_name+' '+ c_index + r_index + '    ' + x + y + z + '  1.00  0.00'+line_feed

        return line


    def build_oniom(self, layer, chrg=None, cnt_info:list=None, if_lig=0 , if_sol=0):
        '''
        build line for oniom. Use element name for ligand atoms since they are mostly in QM regions.
        Gaussian use *ff96* which the atom type is corresponding to all_amino94.lib and ion94.lib in Amber distribution
        ---------
        chrg    : charge of current atom
        layer   : layer where atom in
        cnt_info: (low atom only) when current atom was a boundary atom. Provide [cnt_atom_ele, cnt_atom_label, cnt_atom_id]
        '''
        cnt_flag = ''
        self.get_ele()
        if layer not in ['h','l']:
            raise Exception('build_oniom: please use: "h" or "l" for layer')
        if layer == 'h':
            fz_flag = '0'
            ly_flag = 'H'
        if layer == 'l':
            fz_flag = '-1'
            ly_flag = 'L'       
            if cnt_info != None:
                cnt_flag = ' '+ cnt_info[0] + '-'+cnt_info[1]+' '+str(cnt_info[2])

        # label: deal with N/C terminal and ligand            
        if if_lig:
            G16_label = self.ele
        if if_sol:
            G16_label = G16_label_map[self.parent.name][self.name]
        if not if_sol and not if_lig:
            r = self.resi
            r1 = self.resi.chain.residues[0]
            rm1 = self.resi.chain.residues[-1]
            if r == r1 and self.name in ['H1','H2','H3']:
                G16_label = 'H'
            else:
                if r == rm1 and self.name == 'OXT':
                    G16_label = 'O2'
                else:
                    G16_label = G16_label_map[self.parent.name][self.name]

        # chrg
        if chrg == None:
            try:
                chrg = self.charge # from prmtop
            except NameError:
                raise Exception('You need to at least provide a charge or use get_atom_charge to get one from prmtop file.')

        atom_label = '{:<16}'.format(' '+self.ele+'-'+G16_label+'-'+str(round(chrg,6)))
        fz_flag = '{:>2}'.format(fz_flag)
        x = '{:<14.8f}'.format(self.coord[0])
        y = '{:<14.8f}'.format(self.coord[1])
        z = '{:<14.8f}'.format(self.coord[2])

        line = atom_label+' '+ fz_flag + '   ' + x +' '+ y +' '+ z + ' ' + ly_flag + cnt_flag +line_feed

        return line


    def get_around(self, rad):
        pass

    
    def get_ele(self):
        '''
        get self.ele from a certain map according to the ff type
        '''
        if self.name in Resi_Ele_map[self.ff].keys():
            self.ele = Resi_Ele_map[self.ff][self.name]
        else:
            self.ele = re.match('^[A-Z][a-z]?',self.name).group()

    '''
    ====
    Special Method
    ====
    '''

    def __getattr__(self, key):
        if key == 'residue' or key == 'resi':
            return self.parent
        else:
            Exception('bad key: getattr error')
    
    def __int__(self):
        return self.id



class Metalatom(Atom):
    '''
    -------------
    initilize from
    PDB:        Atom.fromPDB(atom_input, input_type='PDB_line' or 'line_str' or 'file' or 'path')
    raw data:   Atom(atom_name, coord)
    Atom:       MetalAtom.fromAtom(atom_obj)
    -------------
    id
    name
    coord = [x, y, z]
    ele
    resi_name

    parent # the whole stru
    -------------
    Method
    -------------
    set_parent
    get_donor_atom(self, method='INC', check_radius=4.0)
    get_donor_residue(self, method='INC')
    -------------
    '''


    def __init__(self, name, resi_name, coord, ff, id=None, parent=None):
        '''
        Have both atom_name, ele and resi_name 
        '''
        self.resi_name = resi_name
        self.ele = Metal_map[resi_name] 
        Atom.__init__(self, name, coord, ff, id, parent)

        self.donor_atoms = []
        self.parm=None

    @classmethod
    def fromAtom(cls, atom_obj):
        '''
        generate from Atom object. copy data.
        '''
        return cls(atom_obj.name, atom_obj.parent.name, atom_obj.coord, atom_obj.ff, id=atom_obj.id, parent=atom_obj.parent)

    def get_valence(self):
        pass

    # fix related
    def get_donor_atom(self, method='INC', check_radius=4.0):
        '''
        Get coordinated donor atom for a metal center.
        1. check all atoms by type, consider those in the "donor_map"
        (only those atoms from residues are considered. So atoms from ligand and ion will be ignored)
        2. check if atoms are within the check_radius
        3. check distance for every atoms left.
        -----------------------
        Base on a certain type of radius of this metal:
        method = INC ------ {ionic radius} for both metal and donor atom
               = VDW ------ {Van Der Waals radius} for both metal and donor atom
        check_radius ------ the radius that contains all considered potential donor atoms. set for reducing the complexity.
        -----------------------
        save found atom list to self.donor_atoms
        '''

        # find radius for matal
        if method == 'INC':
            R_m = Ionic_radius_map[self.ele]
        if method == 'VDW':
            R_m = VDW_radius_map[self.ele]
        
        # get target with in check_radius (default: 4A)
        coord_m = np.array(self.coord)
        protein_atoms = self.parent.get_all_protein_atom()
        for atom in protein_atoms:
            
            #only check donor atom (by atom name)
            if atom.name in Donor_atom_list[atom.ff]:
                
                # cut by check_radius
                dist = np.linalg.norm(np.array(atom.coord) - coord_m)
                if dist <= check_radius:
                    # determine coordination
                    atom.get_ele()
                    if method == 'INC':
                        R_d = Ionic_radius_map[atom.ele]
                    if method == 'VDW':
                        R_d = VDW_radius_map[atom.ele]
                    
                    if dist <= (R_d + R_m):
                        self.donor_atoms.append(atom)
                        if Config.debug > 1:
                            print('Metalatom.get_donor_atom: '+self.name+' find donor atom:' + atom.resi.name +' '+ str(atom.resi.id) + ' ' + atom.name)                     
        

    def get_donor_residue(self, method='INC'):
        '''
        get donor residue based on donor atoms
        '''
        self.donor_resi =  []
        self.get_donor_atom(method=method)
        # add d_atom and save donor_resi
        for atom in self.donor_atoms:
            resi = atom.resi
            resi.d_atom = atom
            resi.a_metal = self
            self.donor_resi.append(resi)

        # warn if more than one atom are from a same residue
        for index in range(len(self.donor_resi)):
            for index2 in range(len(self.donor_resi)):
                if index2 > index:
                    if self.donor_resi[index2].id == self.donor_resi[index].id:
                        print('\033[1;31;40m!WARNING! found more than 1 donor atom from residue: '+ self.donor_resi[index].name + str(self.donor_resi[index].id) +'\033[m')


    def _metal_fix_1(self):
        '''
        Fix1: deprotonate all donor (rotate those with tight H, like Ser)
            - warn if uncommon donor residue shows up (like Ser)
        '''
        for resi in self.donor_resi:
            if resi.ifDeProton():
                resi.deprotonate(resi.d_atom)
            else:
                if resi.name not in NoProton_list:
                    print('!WARNING!: uncommon donor residue -- '+resi.chain.id+' '+resi.name+str(resi.id))
                    #resi.rot_proton(resi.d_atom)


    def _metal_fix_2(self):
        '''
        Fix2: rotate if there're still lone pair left 
        '''
        for resi in self.donor_resi:
            resi.rot_proton(resi.d_atom)

    def _metal_fix_3(self):
        '''
        Fix3: run pka calculate containing ion (maybe pypka) and run Fix2 based on the result
        waiting for response
        '''
        pass


    def build(self, a_id=None, r_id=None, c_id=None,  ff='AMBER', forcefield = 'ff14SB'):
        '''
        generate an metal atom output line. End with LF
        return a line str
        --------
        use self.id if not assigned
        
        '''
        #default
        if a_id == None:
            a_id = self.id
        if r_id == None:
            r_id = self.resi.id
        if c_id == None:
            print('WARNING: the metal atom may need a chain id in build()!')
            c_id = ' '
        
        if ff == 'AMBER':
            #build an amber style line here
            l_type = '{:<6}'.format('ATOM')
            a_index = '{:>5d}'.format(a_id)

            if forcefield == 'ff14SB':
                # ff14SB by default
                a_name = '{:>2}'.format(self.name)
                r_name = '{:<3}'.format(self.resi_name)
            else:
                raise Exception('Only support ff14SB atom/resiude name now')

            c_index = c_id
            r_index = '{:>4d}'.format(r_id)
            x = '{:>8.3f}'.format(self.coord[0])
            y = '{:>8.3f}'.format(self.coord[1])
            z = '{:>8.3f}'.format(self.coord[2])
        
        #example: ATOM   5350  HB2 PRO   347      32.611  15.301  24.034  1.00  0.00
        line = l_type + a_index +' '+a_name +'   '+ r_name+' '+ c_index + r_index + '    ' + x + y + z + '  1.00  0.00'+line_feed

        return line


    def build_oniom(self, layer, chrg=None, cnt_info:list=None):
        '''
        build a metal line for oniom. 
        Gaussian use *ff96* which the atom type is corresponding to all_amino94.lib and ion94.lib in Amber distribution
        For metals that do not exist in ion94.lib. Custom a type and get parms from TIP3P lib.
        ---------
        chrg    : charge of current metal atom
        layer   : layer where atom in
        cnt_info: (low atom only) when current atom was a boundary atom. Provide [cnt_atom_ele, cnt_atom_label, cnt_atom_id]
        '''
        cnt_flag = ''
        if layer not in ['h','l']:
            raise Exception('build_oniom: please use: "h" or "l" for layer')
        if layer == 'h':
            fz_flag = '0'
            ly_flag = 'H'
        if layer == 'l':
            fz_flag = '-1'
            ly_flag = 'L'       
            if cnt_info != None:
                cnt_flag = ' '+ cnt_info[0] + '-'+cnt_info[1]+' '+str(cnt_info[2])

        # label
        if self.resi_name in G16_label_map.keys():
            G16_label = G16_label_map[self.resi_name][self.name]
        else:
            if Config.debug >= 1:
                print('Metal: '+self.name+ ' not in build-in atom type of ff96.')
                print('Use parameters and atom types from TIP3P (frcmod.ionsjc_tip3p & frcmod.ions234lm_126_tip3p)')
            G16_label = self.resi_name.strip('+-')+'0'
            self.parm = tip3p_metal_map[self.resi_name][self.name]
            self.parm[0] = G16_label

        # chrg
        if chrg == None:
            try:
                chrg = self.charge # from prmtop
            except NameError:
                raise Exception('You need to at least provide a charge or use get_atom_charge to get one from prmtop file.')

        atom_label = '{:<16}'.format(' '+self.ele+'-'+G16_label+'-'+str(round(chrg,6)))
        fz_flag = '{:>2}'.format(fz_flag)
        x = '{:<14.8f}'.format(self.coord[0])
        y = '{:<14.8f}'.format(self.coord[1])
        z = '{:<14.8f}'.format(self.coord[2])

        line = atom_label+' '+ fz_flag + '   ' + x +' '+ y +' '+ z + ' ' + ly_flag + cnt_flag +line_feed

        return line


    '''
    ====
    Special Method
    ====
    '''

    def __getitem__(self, key: int):
        '''
        Metalatom_obj[0]: self
        -----
        return self when iter. Allow metalatom to be iterate like a residue that return a atom level obj.
        '''
        if key == 0:
            return self
        if key > 0:
            raise StopIteration



class Ligand(Residue):
    '''
    -------------
    initilize from
    PDB:        Residue.fromPDB(atom_input, input_type='PDB_line' or 'line_str' or 'file' or 'path')
    raw data:   Residue(atom_name, coord)
    Residue:    Ligand.fromResidue(atom_obj)
    -------------
    id
    name
    atoms = [atom, ...]
    parent # the whole stru
    -------------
    Method
    -------------
    set_parent
    -------------
    '''
    def __init__(self, atoms, id, name, net_charge=None, parent=None):
        self.net_charge = net_charge
        Residue.__init__(self, atoms, id, name, parent)

    @classmethod
    def fromResidue(cls, Resi_obj):
        '''
        generate from a Residue object. copy data
        '''
        return cls(Resi_obj.atoms, Resi_obj.id, Resi_obj.name, parent=Resi_obj.parent)

    @classmethod
    def fromPDB(cls, resi_input, resi_id=None, resi_name=None, net_charge=None, input_type='PDB_line'):
        '''
        generate resi from PDB. Require 'ATOM' and 'HETATM' lines.
        ---------
        resi_input = PDB_line (or line_str or file or path)
        resi_id : int (use the number in the line by default // support customize)
        net_charge : user assigned net charge for further use
        Use PDB_line in the list to init each atom
        '''

        # adapt general input // converge to a list of PDB_line (resi_lines)
        if input_type == 'path':
            f = open(resi_input)
            resi_lines = PDB_line.fromlines(f.read())
            f.close()
        if input_type == 'file':
            resi_lines = PDB_line.fromlines(resi_input.read())
        if input_type == 'line_str':
            resi_lines = PDB_line.fromlines(resi_input)
        if input_type == 'PDB_line':
            resi_lines = resi_input
        
        #clean lines
        for i in range(len(resi_lines)-1,-1,-1):
            if resi_lines[i].line_type != 'ATOM' and resi_lines[i].line_type != 'HETATM':
                if Config.debug > 1:
                    print('Residue.fromPDB: delete error line in input.')
                del resi_lines[i]
    
        # Default resi_id
        if resi_id is None:
            resi_id = resi_lines[0].resi_id
        # get name from first line
        if resi_name is None:
            resi_name = resi_lines[0].resi_name
        # get child atoms
        atoms = []
        for pdb_l in resi_lines:
            atoms.append(Atom.fromPDB(pdb_l))
        
        return cls(atoms, resi_id, resi_name, net_charge=net_charge)
    
 
    def sort(self):
        pass


    def build(self, out_path, ft='PDB'):
        '''
        build ligand file(out_path) with ft format
        '''
        if ft == 'PDB':
            with open(out_path,'w') as of:
                a_id = 0
                for atom in self:
                    a_id = a_id + 1
                    line = atom.build(a_id = a_id, c_id=' ')
                    of.write(line)
                of.write('TER'+line_feed+'END'+line_feed)
        else:
            raise Exception('Support only PDB output now.')


    def get_net_charge(self, method='PYBEL', ph=7.0, o_dir='.'):
        '''
        get net charge for the ligand
        -------
        method   : PYBEL (default) use UNITY_ATOM_ATTR info from openbabel mol2
        '''
        # build file
        mkdir(o_dir+'/cache')
        temp_path=o_dir+'/cache/ligand_temp.pdb'
        temp_pdb2_path=o_dir+'/cache/ligand_temp2.pdb'
        temp_pdb3_path=o_dir+'/cache/ligand_temp3.pdb'
        temp_m2_path=o_dir+'/cache/ligand_temp.mol2'
        self.build(temp_path)
        # charge
        if method == 'PYBEL':
            pybel.ob.obErrorLog.SetOutputLevel(0)
            # remove H (or the )
            mol = next(pybel.readfile('pdb', temp_path))
            mol.removeh()
            mol.write('pdb', temp_pdb2_path, overwrite=True)
            # clean connectivity
            os.system('cat '+temp_pdb2_path+' |grep \'ATOM\' > '+temp_pdb3_path)
            # add H and result net charge
            mol = next(pybel.readfile('pdb', temp_pdb3_path))
            mol.OBMol.AddHydrogens(False, True, ph)
            mol.write('mol2', temp_m2_path, overwrite=True)
            mol = next(pybel.readfile('mol2', temp_m2_path))
            net_charge=0
            for atom in mol:
                net_charge=net_charge+atom.formalcharge
        self.net_charge=net_charge
        return net_charge

        

class Solvent(Residue):
    '''
    -------------
    initilize from
    PDB:        Residue.fromPDB(atom_input, input_type='PDB_line' or 'line_str' or 'file' or 'path')
    raw data:   Residue(atom_name, coord)
    Residue:    Ligand.fromResidue(atom_obj)
    -------------
    id
    name
    atoms = [atom, ...]
    parent # the whole stru
    -------------
    Method
    -------------
    set_parent
    -------------
    '''
    def __init__(self, atoms, id, name, parent=None):
        Residue.__init__(self, atoms, id, name, parent)

    @classmethod
    def fromResidue(cls, Resi_obj):
        '''
        generate from a Residue object. copy data
        '''
        return cls(Resi_obj.atoms, Resi_obj.id, Resi_obj.name, parent=Resi_obj.parent)

    
    def sort(self):
        pass

