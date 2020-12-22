__doc__='''
This module extract and operate structural infomation from PDB
# will replace some local function in PDB class in the future.
-------------------------------------------------------------------------------------
Class Structure
-------------------------------------------------------------------------------------
__init__(self,path)

'''

class structure(object):

    def __init__(self, PDB_obj):
        '''
        extract the structure from PDB_obj.path
        Target:
        - coordinate
        - connectivity
        - [Chain - Residue - Atom] Tree
        - ... (add upon usage)
        ''' 
        self.path=PDB_obj.conf_path
        self._get_structure()
    
    def _get_structure(self):
        '''
        get structure from the PDB file.
        '''
    def build(self, ff='AMBER'):
        '''
        build PDB after the change
        based on atom and resinames
        '''
        pass

    def find_metal(self):
        '''
        find metal atom return a list of MetalAtom object
        '''
        return []




class chain(structure):
    pass
class residue(chain):
    def deprotonate(self):
        pass

class atom(residue):
    def get_line(self, ff='AMBER'):
        pass
    def get_around(self, rad):
        pass
class MetalAtom(atom):
    def get_valance(self):
        pass
    def get_donor_atom(self):
        pass
    def get_donor_residue(self):
        return []
    

