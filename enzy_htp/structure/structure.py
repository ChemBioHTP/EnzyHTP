"""Defines a structure class that serves as a parent class for different types of proteins.

Author: Chris Jurich, <chris.jurich@vanderbilt.edu>
Date: 2022-02-01
"""
import logging
from biopandas.pdb import PandasPdb


class Structure(PandasPdb):

    def __init__(self):
        # TODO make residues and chains as well
        pass


#		self.sequence_ = kwargs.get('sequence', None)
#		self.structure_ = kwargs.get('structure', None)

    def __eq__(self, other):
        return self.sequence_ == other.sequence_

    def rm_wat(self):
        # TODO put a deprecation warning
        self.rm_solvent()

    def rm_solvent(self):
        print(self.df.keys())

        def clean_df(key_name):
            SKIP = ["Na+", "Cl-", "WAT", "HOH"]
            sliced_df = self.df[key_name]
            print(sliced_df)
            self.df[key_name] = sliced_df[(~sliced_df.atom_name.isin(SKIP)) &
                                          (sliced_df.record_name != 'TER')&(~sliced_df.residue_name.isin(["WAT"]))]

        clean_df('ATOM')
        #clean_df( 'OTHERS' )
        temp = self.df['OTHERS']
        keep = temp.record_name != 'CRYST1'
        self.df['OTHERS'] = temp[keep]
    
    def get_protonation(self):
        pass
