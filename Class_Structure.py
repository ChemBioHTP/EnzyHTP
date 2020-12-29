from Class_line import PDB_line
from helper import Child, line_feed
from AmberMaps import *
__doc__='''
This module extract and operate structural infomation from PDB
# will replace some local function in PDB class in the future.
-------------------------------------------------------------------------------------
Class Structure
-------------------------------------------------------------------------------------
__init__(self,path)

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
    ------------

    if_metalatom

    # Move in in the future
    if_art_resi
    if_ligand
    if_complete
    ligand
    '''

    '''
    ====
    init
    ====
    '''
    ifsort = 0

    def __init__(self, chains=None, metalatoms=None, ligands=None):
        '''
        Common part of init methods: direct from data objects
        '''
        if chains is None and metalatoms is None and ligands is None:
            raise ValueError('need at least one input')
            
        self.chains = []
        self.metalatoms = []
        self.ligands = []

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
    
    @classmethod
    def fromPDB(cls, input_obj, input_type='path'):
        '''
        extract the structure from PDB path. Capable with raw experimental and Amber format
        ---------
        input = path (or file or file_str)
        split the file_str to chain and init each chain
        Target:
        - structure - chain - residue - atom
                   |- metalatom(atom)
                   |- ligand(residue)
        - ... (add upon usage)
        ---------
        Special Method
        ---------
        __len__
            len(obj) = len(obj.child_list)
        ''' 

        # adapt general input // converge to file_str
        if input_type == 'path':
            f = open(input_obj)
            file_str = f.read()
            f.close
        if input_type == 'file':
            file_str = input_obj.read()
        if input_type == 'file_str':
            file_str = input_obj
        
        raw_chains = []
        # get raw chains
        chains_str = file_str.split(line_feed+'TER') # Note LF is required
        for index, chain_str in enumerate(chains_str):
            if chain_str.strip() != 'END':
                Chain_index = chr(65+index) # Covert to ABC using ACSII mapping
                # Generate chains
                raw_chains.append(Chain.fromPDB(chain_str, Chain_index))
        
        # clean chains
        # clean metals
        raw_chains_woM, metalatoms = cls._get_metalatoms(raw_chains, method='1')
        # clean ligands
        raw_chains_woM_woL, ligands = cls._get_ligand(raw_chains_woM)

        return cls(raw_chains_woM_woL, metalatoms, ligands)


    @classmethod
    def _get_metalatoms(cls, raw_chains, method='1'):
        '''
        get metal from raw chains and clean chains by deleting the metal part
        -----
        Method 1:   Assume metal/ligand/solvent can be in any chain
                    Assume all resiude have unique index. 
        (Slow but general)
        
        Method 2: Assume metal/ligand/solvent can only be in a seperate chain
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
                        del chain[residue]
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
    def _get_ligand(cls, raw_chains):
        '''
        get ligand from self.chains and clean chains by deleting the ligand part
        -----
        (TO BE DONE)
        '''
        ligands = []
        # get ligand list
        # convert to ligand obj
        # clean empty chains
        return raw_chains, ligands

    '''
    ====
    Methods
    ====
    '''
    def get_metal_center(self):
        '''
        Extract metal centers from metalatoms. Judged by the MetalCenter_map
        save to self.metal_centers
        return a bool: if there's metal_center
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
        '''
        # list
        if type(obj) == list:
            
            obj_ele=obj[0]

            if type(obj_ele) != Chain and type(obj_ele) != Metalatom and type(obj_ele) != Ligand:
                raise TypeError('Add() method only take Chain / Metalatom / Ligand')

            # add parent and clear id (if sort) assign id (if assigned) leave mark if sort and assigned
            for i in obj:               
                i.Add_parent(self)
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

        # single building block
        else:
            if type(obj) != Chain and type(obj) != Metalatom and type(obj) != Ligand:
                raise TypeError('Add() method only take Chain / Metalatom / Ligand')
            
            obj.Add_parent(self)
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
            
        if sort:
            self.sort()
            

        

    def sort(self):
        '''
        assign index according to current items
        chain.id
        resi.id
        atom.id
        '''
        self.ifsort = 1
        pass

    def build(self, path, ff='AMBER'):
        '''
        build PDB after the change
        based on atom and resinames
        '''
        pass

    def protonation_metal_fix(self, Fix):
        pass

    '''
    ====
    Special Method
    ====
    '''

    def __len__(self):
        '''
        len(obj) = len(obj.child_list)
        '''
        return len(self.chains)+len(self.metalatoms)+len(self.ligands)





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
    Add_parent
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
    self.ifsort = 0
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
            chain_str = f.read()
            f.close()
        if input_type == 'file':
            chain_str = chain_input.read()
        if input_type == 'file_str':
            chain_str = chain_input


        # chain residues
        residues = []
        resi_lines = [] # data holder
        lines = PDB_line.fromlines(chain_str) # Note LF is required
        for i, pdb_l in enumerate(lines):
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
                
                # Deal with the last residue
                if i == len(lines)-1:
                    last_resi = Residue.fromPDB(resi_lines, pdb_l.resi_id)
                    residues.append(last_resi)
                
                # Update for next loop                
                last_resi_index = pdb_l.resi_id

        return cls(residues, chain_id)

    '''
    ====
    Method
    ====
    '''
    def add(self, obj, id=None):
        '''
        1. judge obj type
        2. clean original id
        3. add to corresponding list
        '''
        pass

    def sort(self):
        '''
        maybe useful with other format
        ----
        assign index according to current items
        resi.id
        atom.id
        '''
        self.ifsort = 1
        pass


    def get_chain_seq(self):
        pass

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
        Chain_obj.123 = Chain_obj.residues[123-1] // index mimic (start from 1)
        Chain_obj.HIS = Chain_obj.find_resi_name('HIS') // search mimic
        '''
        if type(key) == int:
            return self.residues[key-1]
        if type(key) == str:
            return self._find_resi_name(key)
        if key == 'stru':
            return self.parent
        Exception('bad key: getattr error')


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
    
    #TODO
    if_art_resi
    -------------
    Method
    -------------
    Add_parent
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
    def deprotonate(self):
        pass

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

    def add(self, obj, id=None):
        '''
        1. judge obj type
        2. clean original id
        3. add to corresponding list
        '''
        pass

    def sort(self):
        '''
        None
        '''
        pass



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
        Residue_obj.123 = Residue_obj.atoms[123-1] // index mimic (start from 1)
        Residue_obj.CA = Residue_obj.find_atom_name('CA') // search mimic
        '''
        if type(key) == int:
            return self.atoms[key-1]
        if type(key) == str:
            return self._find_atom_name(key)
        if key == 'chain':
            return self.parent
        Exception('bad key: getattr error')

    
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

    parent # the whole chain
    -------------
    Method
    -------------
    Add_parent
    Add_id # use only after the whole structure is constructed
    -------------
    '''

    def __init__(self, atom_name: str, coord: list, atom_id = None, parent = None):
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
        self.atom_id = atom_id
    
    @classmethod
    def fromPDB(cls, atom_input, input_type='PDB_line'):
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

        return cls(atom_name, coord)


    '''
    ====
    Method
    ====
    '''

    def gen_line(self, ff='AMBER'):
        '''
        generate an output line. End with LF
        -------
        must use after sort!
        '''
        a_index = '{:>5d}'.format(self.id)
        a_name = '{:<4}'.format(self.name)
        r_name = '{:>3}'.format(self.parent.name) # fix for metal
        line = 'ATOM  '+ a_index +' '+a_name+r_name'    348      21.367   9.559  -3.548  1.00  0.00'+line_feed
        待办
        return line

    def get_around(self, rad):
        pass


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
    

class Metalatom(Atom):

    def __init__(self, name, resi_name, coord, id=None, parent=None):
        '''
        Have both atom_name and resi_name
        '''
        self.resi_name = resi_name
        Atom.__init__(self, name, coord, id, parent)

    @classmethod
    def fromAtom(cls, atom_obj):
        '''
        generate from Atom object. copy data.
        '''
        return cls(atom_obj.atom_name, atom_obj.parent.name, atom_obj.coord, parent=atom_obj.parent)

    def get_valance(self):
        pass
    def get_donor_atom(self):
        pass
    def get_donor_residue(self):
        return []

class Ligand(Residue):
    pass
    

