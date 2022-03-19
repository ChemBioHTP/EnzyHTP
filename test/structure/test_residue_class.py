"""Testing the enzy_htp.structure.Residue() class.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
import os
import pandas as pd
from typing import List
from copy import deepcopy
from collections import defaultdict
from biopandas.pdb import PandasPdb

from enzy_htp.structure import Atom, Residue
CURR_DIR = os.path.dirname(os.path.abspath(__file__))

TEST_PDB_FILE = f"{CURR_DIR}/data/3NIR.pdb"

def make_residues( pdbname : str ) -> List[Residue]:
    """Helper method that retrieves a list of residues from a PDB file."""
    holder = defaultdict(list) 
    p_df : pd.DataFrame = PandasPdb().read_pdb( pdbname ).df['ATOM']
    
    for i, row in  p_df.iterrows():
        atom : Atom = Atom(**row)
        holder[atom.residue_key()].append( atom )
    
    residues = []     
    for rkey, atoms in holder.items():
        residues.append( Residue( rkey, atoms ))
    return residues

RESIDUES : List[Residue] = make_residues( TEST_PDB_FILE )

def test_residue_line_functions():
    """Checking that a number of line-related functions for the Residue() class work."""
    for idx, res in enumerate( RESIDUES[:-1] ):
        one_ahead = RESIDUES[idx+1]
        assert res.max_line() + 1 == one_ahead.min_line()
        assert res.line_range()[-1] + 1 == one_ahead.line_range()[0]
        assert res.neighbors( one_ahead )
        assert one_ahead.neighbors( res )

      
    for idx1, res1 in enumerate( RESIDUES ):
        count = 0
        for idx2, res2 in enumerate( RESIDUES ):
            if idx1 == idx2:
                continue
            count += int(res1.neighbors(res2))

        if not idx1 or idx1 == len(RESIDUES)-1:
            assert count == 1 
        else:
            assert count == 2

def test_residue_key_information():
    """Ensuring that general getting and setting of residue key information works."""
    local_res = deepcopy( RESIDUES[0] )
    assert local_res.residue_key == "A.THR.1"
    local_res.set_chain( "B" )
    assert local_res.residue_key == "B.THR.1"
    assert local_res.num() == 1 

    for aa in local_res.atom_list():
        assert aa.chain_id == "B"

def test_check_all_canonical():
    """Checking that all the loaded Residue()'s are canonical amino acids and not metals, rd_solvents of rd_non_ligands"""
    for rr in RESIDUES:
        assert rr.is_canonical()
        assert not rr.is_metal()
        assert not rr.is_rd_solvent()
        assert not rr.is_rd_non_ligand()

#TODO(CJ) add in some tests for other types of residues
