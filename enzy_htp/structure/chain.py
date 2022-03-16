# TODO what the heck else does this class do???
# CJ: I think this class should be able to remove bad residues on its own
# Should this class know what PDB it came from?
from ..core import _LOGGER
from typing import List

from .residue import Residue

class Chain:
    def __init__(self, name, residues):
        self.name_ = name
        self.residues_ : List[Residue] = residues
        # CJ: not sure if I should do this but I am overwriting every time
        # right now
        _ = list(map(lambda r: r.set_chain(self.name_), self.residues_))

    def is_metal(self) -> bool:
        for rr in self.residues_:
            if rr.is_metal():
                _LOGGER.warn(f"Structure: found metal in raw: {self.name_} {rr.name} {residue.id}")
                return True
        return False

    def is_HET(self) -> bool:
        for rr in self.residues_:
            if not rr.is_canonical():
                return False
        return True

    def empty(self) -> bool:
        return len(self.residues_) == 0

    def residues(self) -> List[Residue]:
        return self.residues_
   
    def name(self) -> str:
       return self.name_ 
    def __getitem__(self, key: int) -> Residue:
       return self.residues[key]


    def __delitem__(self, key: int) -> None:
       del self.residues_[key]



#class Chain(Child):
#    """
#    -------------
#    initilize from
#    PDB:        Chain.fromPDB(chain_input, chain_id, input_type='file_str' or 'file' or 'path')
#    raw data:   Chain(residues, chain_id)
#    -------------
#    id
#    parent # the whole structure
#    residues = [resi_obj, ...]
#    chain_seq = ['resi_name', ..., 'NAN', ..., 'resi_name']
#    -------------
#    method
#    -------------
#    set_parent
#    get_chain_seq(self)
#    _find_resi_name
#    -------------
#    Special method
#    -------------
#    __getitem__
#        Chain_obj[int]: Chain_obj.residues[int] // (start from 0)
#    __getattr__
#        Chain_obj.123 = Chain_obj.residues[123-1] // index mimic (start from 1)
#        Chain_obj.HIS = Chain_obj.find_resi_name('HIS') // search mimic
#    __delitem__
#        del obj[int] --> obj.child[int].remove() // delete by index (start from 0)
#        del obj[str] --> obj._del_child_name() // delete by value
#        del obj[child] --> obj.child_list.remove(child) // delete by value
#    __len__
#        len(obj) = len(obj.child_list)
#    """
#
#    """
#    ====
#    init
#    ====
#    """
#
#    def __init__(self, residues, chain_id: str, parent=None):
#        """
#        Common part of init methods: direct from data objects
#        No parent by default. Add parent by action
#        """
#        # set parent to None
#        Child.__init__(self)
#        # add parent if provided
#        if parent != None:
#            self.set_parent(parent)
#        # adapt some children
#        self.residues = []
#        for i in residues:
#            i.set_parent(self)
#            self.residues.append(i)
#        # set id
#        self.id = chain_id
#
#        # init
#        self.ifsorted = 0
#
#    @classmethod
#    def fromPDB(cls, chain_input, chain_id, input_type="file_str"):
#        """
#        generate chain from PDB. Capable with raw experimental and Amber format. Only read 'ATOM' and 'HETATM' lines.
#        ---------
#        chain_input = file_str (or path or file)
#        chain_id : str
#        split the file_str to residues and init each residue
#        """
#
#        # adapt general input // converge to file_str
#        if input_type == "path":
#            f = open(chain_input)
#            chain_str = f.read().strip()
#            f.close()
#        if input_type == "file":
#            chain_str = chain_input.read().strip()
#        if input_type == "file_str":
#            chain_str = chain_input.strip()
#
#        # chain residues
#        residues = []
#        resi_lines = []  # data holder
#        lines = PDB_line.fromlines(chain_str)  # Note LF is required
#        for pdb_l in lines:
#            if pdb_l.line_type == "ATOM" or pdb_l.line_type == "HETATM":
#                # Deal with the first residue
#                if len(residues) == 0 and len(resi_lines) == 0:
#                    resi_lines.append(pdb_l)
#                    last_resi_index = pdb_l.resi_id
#                    continue
#
#                # find the beginning of a new residue
#                if pdb_l.resi_id != last_resi_index:
#                    # Store last resi
#                    last_resi = Residue.fromPDB(resi_lines, last_resi_index)
#                    residues.append(last_resi)
#                    # empty the holder for current resi
#                    resi_lines = []
#
#                resi_lines.append(pdb_l)
#
#                # Update for next loop
#                last_resi_index = pdb_l.resi_id
#
#        # Deal with the last residue (Now add the last one after the loop.)
#        # if i == len(lines)-1:
#        if resi_lines != []:  # deal with blank chain
#            last_resi = Residue.fromPDB(resi_lines, resi_lines[-1].resi_id)
#            residues.append(last_resi)
#        else:
#            if Config.debug >= 1:
#                print("Chain.fromPDB: find a empty chain: " + chain_id)
#
#        return cls(residues, chain_id)
#
#    """
#    ====
#    Method
#    ====
#    """
#
#    def add(self, obj, id=None, sort=0):
#        """
#        1. judge obj type
#        2. clean original id
#        3. add to corresponding list
#        """
#        # list
#        if type(obj) == list:
#
#            obj_ele = obj[0]
#
#            if type(obj_ele) != Residue:
#                raise TypeError("chain.Add() method only take Residue")
#
#            # add parent and clean id (if sort) assign id (if assigned) leave mark if sort and assigned
#            for i in obj:
#                i.set_parent(self)
#                if sort:
#                    if id != None:
#                        i.id = str(id) + "i"  # str mark
#                    else:
#                        i.id = id  # None
#                else:
#                    if id != None:
#                        i.id = id
#            self.residues.extend(obj)
#
#        # single building block
#        else:
#            if type(obj) != Residue:
#                raise TypeError("Chain.Add() method only take Residue")
#
#            obj.set_parent(self)
#            if sort:
#                if id != None:
#                    obj.id = str(id) + "i"  # str mark
#                else:
#                    obj.id = id  # None
#            else:
#                if id != None:
#                    obj.id = id
#            self.residues.append(obj)
#
#        if sort:
#            self.sort()
#
#    def sort(self, sort_resi=1):
#        """
#        sort_resi: 1/0 -- if or not sort atoms in each residue. (start from 1 for each residue)
#
#        turn residue index into str and sort with list.sort() -- cope with the insert mark
#        * start form 1 in each chain
#        * if added object has same id and is not assigned with a insert mark -- place after a original one.
#        ----
#        assign index according to current items
#        resi.id
#        atom.id
#        """
#        # sort residue order
#        for i in self.residues:
#            if type(i.id) == str:
#                raise Exception(
#                    "Does not support added residue now: update in the future"
#                )
#            # TODO
#
#        self.residues.sort(key=lambda i: i.id)
#        # re-id each residue
#        for index, resi in enumerate(self.residues):
#            resi.id = index + 1
#            # sort each residue
#            if sort_resi:
#                resi.sort()
#
#        self.ifsorted = 1
#
#    def get_chain_seq(self, Oneletter=0):
#        """
#        get chain sequence from the "residues" list. Store in self.seq
#        ----------------------------
#        Oneletter
#        = 0 // use 3-letter format to represet each residue
#        = 1 // use 1-letter format to represet each residue
#        ----------------------------
#        + Use "NAN"/"-" as a filler to store missing residues (detect internal missing)
#        - (WARNING) Require the ligand in seperate chains
#        - (WARNING) Only support normal chain! (not ligand or solvent)
#        ---------------------
#        * Note that the missing sequence infomation will be missing after sort(). So get seq before sort to obtain such info.
#        * Need to be update after mutation
#        """
#        if self.ifsorted:
#            if Config.debug >= 1:
#                print(
#                    "Chain.get_chain_seq: WARNING: missing sequence infomation is missing after sort()."
#                )
#
#        chain_seq = []
#        for resi in self.residues:
#            # first resi
#            if len(chain_seq) == 0:
#                last_id = resi.id
#                chain_seq.append(resi.name)
#                continue
#            # missing sequence
#            missing_length = resi.id - last_id - 1
#            if missing_length > 0:
#                chain_seq = chain_seq + ["NAN",] * missing_length
#
#            chain_seq.append(resi.name)
#            last_id = resi.id
#
#        self.seq = chain_seq
#        if Oneletter == 1:
#            self._get_Oneletter()
#            return self.seq_one
#        else:
#            return self.seq
#
#    def if_art_resi(self):
#        """
#        TODO
#        """
#        pass
#
#    def _get_Oneletter(self):
#        """
#        (Used internally) convert sequences in self.sequence to oneletter-based str
#        - The 'NAN' is convert to '-'
#        - Present unnature residue as full 3-letter name
#        save to self.seq_one
#        """
#        if len(self.seq) == 0:
#            raise IndexError("The self.sequence should be obtained first")
#
#        seq_one = ""
#        for name in self.seq:
#            if name == "NAN":
#                seq_one = seq_one + "-"
#            else:
#                if name in Resi_map2:
#                    seq_one = seq_one + Resi_map2[name]
#                else:
#                    seq_one = seq_one + " " + name + " "
#
#        self.seq_one = seq_one
#
#    def _find_resi_name(self, name: str):
#        """
#        find residues according to the name
#        return a list of found residues
#        """
#        out_list = []
#        for resi in self.residues:
#            if resi.name == name:
#                out_list.append(resi)
#        return out_list
#
#    def _find_resi_id(self, id: int):
#        """
#        find residues according to the id
#        return a found residue
#        """
#        out_list = []
#        for resi in self.residues:
#            if resi.id == id:
#                out_list.append(resi)
#        if len(out_list) > 1:
#            print(
#                "\033[32;0mShould there be same residue id in chain +"
#                + self.name
#                + str(self.id)
#                + "?+\033[0m"
#            )
#            raise Exception
#        return out_list[0]
#
#    def _del_resi_name(self, name: str):
#        """
#        find residues according to the name
#        delete found residues
#        """
#        for i in range(len(self.residues) - 1, -1, -1):
#            if self.residues[i].name == name:
#                del self.residues[i]
#
#    """
#    ====
#    Special Method
#    ====
#    """
#
#    def __getitem__(self, key: int):
#        """
#        Chain_obj[int]: Chain_obj.residues[int]
#        -----
#        use residue index within the chain, start from 0
#        """
#        return self.residues[key]
#
#    def __getattr__(self, key):
#        """
#        Use i to indicate id (! by this, do not support searching any residue name start with i)
#        Chain_obj.i123 = Chain_obj.residues[123-1] // index mimic (start from 1)
#        Chain_obj.HIS = Chain_obj.find_resi_name('HIS') // search mimic
#        """
#        if key == "stru":
#            return self.parent
#        if key[0] == "i":
#            # digit
#            key = int(key[1:])
#            return self.residues[key - 1]
#        else:
#            # text
#            return self._find_resi_name(key)
#
#    def __delitem__(self, key):
#        """
#        del obj[int] --> obj.child[int].remove() // delete by index (start from 0)
#        del obj[str] --> obj._del_child_name() // delete by value
#        del obj[child] --> obj.child_list.remove(child) // delete by value
#        """
#        if type(key) == int:
#            del self.residues[key]
#        if type(key) == str:
#            self._del_resi_name(key)
#        if type(key) == Residue:
#            self.residues.remove(key)
#
#    def __len__(self):
#        """
#        len(obj) = len(obj.child_list)
#        """
#        return len(self.residues)
#
#

