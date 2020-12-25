from Class_line import PDB_line
from helper import Child
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
            self.chains.append(chain.add_parent(self))
        for metalatom in metalatoms:
            self.metalatoms.append(metalatom.add_parent(self))
        for ligand in ligands:
            self.ligands.append(ligand.add_parent(self))
    
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
        chains_str = file_str.split('\nTER') # Note LF is required
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
        Method 1: Assume metal/ligand/solvent can be in any chain
        (Slow but general)
        
        Method 2: Assume metal/ligand/solvent can only be in a seperate chain
        (fast but limited)
        '''
        if method == '1':
            for chain in raw_chains:
                for residue in chain:
                    for atom in residues:
                        pass #写完getitem之后完成这部分 定义delitem通过索引删除元素
                    
        if method == '2':
            pass

        return raw_chains, metalatoms
        # 做这部分！！
        ###############################################
        


    @classmethod
    def _get_ligand(cls, raw_chains):
        '''
        get ligand from self.chains and clean chains by deleting the ligand part
        -----
        (TO BE DONE)
        '''
        return raw_chains, []

    '''
    ====
    Methods
    ====
    '''

    def get_art_resi(self):
        '''
        find art_resi
        '''
        pass
    
    def add_chain(self, chain_obj):
        '''
        1. set parent
        2. add chain
        '''
        pass

    def build(self, ff='AMBER'):
        '''
        build PDB after the change
        based on atom and resinames
        '''
        pass




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
    
    '''

    '''
    ====
    init
    ====
    '''
    def __init__(self, residues, chain_id: str):
        '''
        Common part of init methods: direct from data objects
        No parent by default. Add parent by action
        '''
        #set parent to None
        Child.__init__(self) 
        #adapt some children
        self.residues = []
        for i in residues:
            self.residues.append(i.add_parent(self))
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
                    residues.append(Residue.fromPDB(resi_lines, last_resi_index))
                    # empty the holder for current resi
                    resi_lines = []

                resi_lines.append(pdb_l)
                # Update for next loop                
                last_resi_index = pdb_l.resi_id

        return cls(residues, chain_id)

    '''
    ====
    Method
    ====
    '''

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
            return self.find_resi_name(key)



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
    '''

    '''
    ====
    init
    ====
    '''
    def __init__(self, atoms, resi_id, resi_name):
        '''
        Common part of init methods: direct from data objects
        No parent by default. Add parent by action
        '''
        #set parent to None
        Child.__init__(self) 
        #adapt some children
        self.atoms = []
        for i in atoms:
            self.atoms.append(i.add_parent(self))
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
            return self.find_atom_name(key)


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

    def __init__(self, atom_name: str, coord: list, atom_id = None):
        '''
        Common part of init methods: direct from data objects
        No parent by default. Add parent by action
        '''
        #set parent to None
        Child.__init__(self)
        # get data
        self.name = atom_name
        self.coord = coord
        # self.atom_id = atom_id
    
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
        generate an output line
        '''
        pass

    def get_around(self, rad):
        pass
    

class Metalatom(Atom):

    @classmethod
    def fromAtom(cls, atom_obj):
        '''
        generate from Atom object. cope data
        '''
        return cls(atom_obj.atom_name, atom_obj.coord)

    def get_valance(self):
        pass
    def get_donor_atom(self):
        pass
    def get_donor_residue(self):
        return []
    

