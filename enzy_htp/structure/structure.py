"""Defines a structure class that serves as a parent class for different types of proteins.

Author: Chris Jurich, <chris.jurich@vanderbilt.edu>
Date: 2022-02-01
"""
import logging
from biopandas.pdb import PandasPdb
from pdb2pqr.main import main_driver as run_pdb2pqr
from pdb2pqr.main import build_main_parser as build_pdb2pqr_parser
from .atom import Atom
from .chain import Chain
from .residue import Residue
from collections import defaultdict


class Structure:
    """Represents a protein. Used for """

    def __init__(self, **kwargs):
        self.chains_ = dict()
        self.residues_ = dict()
        self.pandas_pdb = PandasPdb()
        self.current_pdb = str()
        self.temp_dir = kwargs.get('temp_dir',)

    def rm_wat(self) -> None:
        """Removes water and ions from the structure. Warning: Deprecated. Use Structure.rm_solvent() instead"""
        logging.warning(
            'Structure.rm_wat() is deprecated. Use Structure.rm_solvent()')
        self.rm_solvent()

    def read_pdb(self, pdbFile) -> None:
        #TODO add some kind of file validation here
        self.current_pdb = pdbFile
        self.pandas_pdb.read_pdb(pdbFile)

    def to_pdb(self, pdb_name):
        pass

    def rm_solvent(self):
        """Removes water and ions from the structure as well as CRYST1 lines."""
        print(self.df.keys())

        # TODO fix this
        def clean_df(key_name):
            SKIP = ["Na+", "Cl-", "WAT", "HOH"]
            sliced_df = self.df[key_name]
            print(sliced_df)
            self.df[key_name] = sliced_df[
                (~sliced_df.atom_name.isin(SKIP)) &
                (sliced_df.record_name != 'TER') &
                (~sliced_df.residue_name.isin(["WAT"]))]

        clean_df('ATOM')
        temp = self.df['OTHERS']
        keep = temp.record_name != 'CRYST1'
        self.df['OTHERS'] = temp[keep]

    def get_protonation(self, ph=7.0):
        #TODO put in a check for the pH range
        #TODO should probably save this stuff
        #TODO should check that the ff is correct
		#TODO should allow user to specify the force field
        pdb2pqr_parser = build_pdb2pqr_parser()
        load_path = 'output.pdb'
        self.to_pdb(load_path)
        args = pdb2pqr_parser.parse_args([
            "--ff=PARSE",
            "--ffout=" + "AMBER",
            "--with-ph=" + str(ph),
            load_path,
            'output.pqr',
        ])
        run_pdb2pqr(args)
		# basically the idea here is that we want to set the pdb2pqr structure
		# as the new structure, but I guess that sometimes it will drop the ligands

    def build_chains(self):
        """"""
        self.build_residues()
        mapper = defaultdict(list)
        for res in self.residues_.values():
            mapper[res.chain()].append(res)

        for chain_name, residues in mapper.items():
            self.chains_[chain_name] = Chain(chain_name, sorted(residues, key=lambda r: r.num()))

    def build_residues(self) -> None:
        mapper = defaultdict(list)
        for i, row in self.pandas_pdb.df['ATOM'].iterrows():
            aa = Atom(**row)
            mapper[aa.residue_key()].append(aa)

        for res_key, atoms in mapper.items():
            self.residues_[res_key] = Residue(residue_key=res_key, atoms=sorted(atoms, key=lambda a: a.atom_number))
