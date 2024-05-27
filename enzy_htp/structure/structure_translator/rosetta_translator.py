"""Instantiation of TranslatorBase for the Rosetta molecular modelling suite. This class 
enables support for translation of canonical amino acid names to and from the standard
AmberMD scheme in EnzyHTP.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2024-05-27
"""

from plum import dispatch

from ..atom import Atom
from ..residue import Residue
from ..structure import Structure


from .translator_base import TranslatorBase


class RosettaTranslator(TranslatorBase):
    """Instantiation of TranslatorBase() for the Rosetta molecular modellig suite. """

    def naming_style(self) -> str:
        """This is the Rosetta naming scheme."""
        return "rosetta"


    @dispatch
    def to_standard(self, res:Residue) -> None:
        """Special overloaded of to_standard() that addresses funky histidine 3 letter naming."""

        if res.name == 'HIS':
            atom_names = [aa.name for aa in res.atoms]

            has_delta:bool = 'HD1' in atom_names
            has_eps:bool = 'HE2' in atom_names

            new_name:str = None

            if has_delta and has_eps:
                new_name = 'HIP'
            elif has_delta:
                new_name = 'HID'
            elif has_eps:
                new_name = 'HIE'
            else:
                raise ValueError()
            res.name = new_name
        else:
            super().to_standard(res)    


    def init_mappings(self) -> None:
        """Registers mappings for Rosetta to/from AmberMD. Most changes relate to Hydrogens."""
        self.register_mapping('ALA', ['HB1', 'HB2', 'HB3', 'H1', 'H2', 'H3'], 
                                'ALA', ['1HB', '2HB', '3HB', '1H', '2H', '3H'])
        self.register_mapping('ARG', ['HB3', 'HB2', 'HG3', 'HG2', 'HD3', 'HD2', 'HH11', 'HH12', 'HH21', 'HH22', 'H1', 'H2', 'H3'],
                                'ARG', ['1HB', '2HB', '1HG', '2HG', '1HD', '2HD', '1HH1', '2HH1', '1HH2', '2HH2', '1H', '2H', '3H'])
        self.register_mapping('ASN', ['HB3', 'HB2', 'HD21', 'HD22', 'H1', 'H2', 'H3'],
                                'ASN', ['1HB', '2HB', '1HD2', '2HD2', '1H', '2H', '3H'])
        self.register_mapping('ASP', ['HB3', 'HB2', 'H1', 'H2', 'H3'], 
                                'ASP', ['1HB', '2HB', '1H', '2H', '3H'])
        self.register_mapping('CYS', ['HB3', 'HB2', 'H1', 'H2', 'H3'],
                                'CYS', ['1HB', '2HB', '1H', '2H', '3H'])
        self.register_mapping('GLN', ['HB3', 'HB2', 'HG3', 'HG2', 'HE21', 'HE22', 'H1', 'H2', 'H3'], 
                                'GLN', ['1HB', '2HB', '1HG', '2HG', '1HE2', '2HE2', '1H', '2H', '3H'])
        self.register_mapping('GLU', ['HB3', 'HB2', 'HG3', 'HG2', 'H1', 'H2', 'H3'], 
                                'GLU', ['1HB', '2HB', '1HG', '2HG', '1H', '2H', '3H'])
        self.register_mapping('GLY', ['HA3', 'HA2', 'H1', 'H2', 'H3'], 
                                'GLY', ['1HA', '2HA', '1H', '2H', '3H'])
      
        self.register_mapping('HID', ['HB3', 'HB2', 'H1', 'H2', 'H3'], 
                                'HIS', ['1HB', '2HB', '1H', '2H', '3H'])
        self.register_mapping('HIE', ['HB3', 'HB2', 'H1', 'H2', 'H3'], 
                                'HIS', ['1HB', '2HB', '1H', '2H', '3H'])
        self.register_mapping('HIP', ['HB3', 'HB2', 'H1', 'H2', 'H3'], 
                                'HIS', ['1HB', '2HB', '1H', '2H', '3H'])

        self.register_mapping('ILE', ['HG13', 'HG12', 'HG21', 'HG22', 'HG23', 'HD11', 'HD12', 'HD13', 'H1', 'H2', 'H3'], 
                                'ILE', ['1HG1', '2HG1', '1HG2', '2HG2', '3HG2', '1HD1', '2HD1', '3HD1', '1H', '2H', '3H'])
        self.register_mapping('LEU', ['HB3', 'HB2', 'HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23', 'H1', 'H2', 'H3'], 
                                'LEU', ['1HB', '2HB', '1HD1', '2HD1', '3HD1', '1HD2', '2HD2', '3HD2', '1H', '2H', '3H'])
        self.register_mapping('LYS', ['HB3', 'HB2', 'HG3', 'HG2', 'HD3', 'HD2', 'HE3', 'HE2', 'HZ1', 'HZ2', 'HZ3', 'H1', 'H2', 'H3'], 
                                'LYS', ['1HB', '2HB', '1HG', '2HG', '1HD', '2HD', '1HE', '2HE', '1HZ', '2HZ', '3HZ', '1H', '2H', '3H'])
        self.register_mapping('MET', ['HB3', 'HB2', 'HG3', 'HG2', 'HE1', 'HE2', 'HE3', 'H1', 'H2', 'H3'], 
                                'MET', ['1HB', '2HB', '1HG', '2HG', '1HE', '2HE', '3HE', '1H', '2H', '3H'])
        self.register_mapping('PHE', ['HB3', 'HB2', 'H1', 'H2', 'H3'], 
                                'PHE', ['1HB', '2HB', '1H', '2H', '3H'])
        self.register_mapping('PRO', ['HB3', 'HB2', 'HG3', 'HG2', 'HD3', 'HD2', 'H1', 'H2', 'H3'], 
                                'PRO', ['1HB', '2HB', '1HG', '2HG', '1HD', '2HD', '1H', '2H', '3H'])
        self.register_mapping('SER', ['HB3', 'HB2', 'H1', 'H2', 'H3'], 
                                'SER', ['1HB', '2HB', '1H', '2H', '3H'])
        self.register_mapping('THR', ['HG21', 'HG22', 'HG23', 'H1', 'H2', 'H3'], 
                                'THR', ['1HG2', '2HG2', '3HG2', '1H', '2H', '3H'])
        self.register_mapping('TRP', ['HB3', 'HB2', 'H1', 'H2', 'H3'], 
                                'TRP', ['1HB', '2HB', '1H', '2H', '3H'])
        self.register_mapping('TYR', ['HB3', 'HB2', 'H1', 'H2', 'H3'], 
                                'TYR', ['1HB', '2HB', '1H', '2H', '3H'])
        self.register_mapping('VAL', ['HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23', 'H1', 'H2', 'H3'], 
                                'VAL', ['1HG1', '2HG1', '3HG1', '1HG2', '2HG2', '3HG2', '1H', '2H', '3H'])

