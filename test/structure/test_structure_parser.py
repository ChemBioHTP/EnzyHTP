"""Testing the functions in the enzy_htp.structure.structure_parser.py file. These functions load and distinguish different types of structural
elements in the enzy_htp framework.

Author: Chris Jurich <chris.jurich@vanderbilt.edu
Date: 2022-03-31
"""
import os
import pytest
import string
import warnings
import pandas as pd
from typing import Dict
from biopandas.pdb import PandasPdb

import enzy_htp.chemical as chem
from enzy_htp.core import file_system as fs
from enzy_htp.structure import structure_parser as sp
from enzy_htp.structure import (
    Structure,
    Residue,
    Chain,
    structure_from_pdb,
    Atom,
    MetalAtom,
    Ligand,
    Solvent,
)


CURRDIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = f"{CURRDIR}/data/"


def test_check_valid_pdb_good_input():
    """Good input for the check_valid_pdb() helper method."""
    dummy_pdb = f"{CURRDIR}/dummy.pdb"
    assert not os.path.exists(dummy_pdb)
    fs.write_lines(dummy_pdb, ["line1", "line2"])
    assert not sp.check_valid_pdb(dummy_pdb)
    fs.safe_rm(dummy_pdb)
    assert not os.path.exists(dummy_pdb)


def test_check_valid_pdb_bad_input():
    """Good input for the check_valid_pdb() helper method."""
    # non pdb file
    txt_file = "not_pdb.txt"
    assert not os.path.exists(txt_file)
    with pytest.raises(SystemExit) as exe:
        sp.check_valid_pdb(txt_file)

    assert exe
    assert exe.type == SystemExit
    assert exe.value.code == 1

    # pdb that doesn't exist
    pdb_imaginary = "not_real.pdb"
    with pytest.raises(SystemExit) as exe:
        sp.check_valid_pdb(pdb_imaginary)

    assert exe
    assert exe.type == SystemExit
    assert exe.value.code == 1

    non_ascii_pdb = f"{CURRDIR}/bad_pdb.pdb"
    assert not os.path.exists(non_ascii_pdb)
    fs.write_lines(non_ascii_pdb, ["日本人 中國的"])
    assert os.path.exists(non_ascii_pdb)

    with pytest.raises(SystemExit) as exe:
        sp.check_valid_pdb(non_ascii_pdb)

    assert exe
    assert exe.type == SystemExit
    assert exe.value.code == 1

    fs.safe_rm(non_ascii_pdb)
    assert not os.path.exists(non_ascii_pdb)


def test_legal_chain_names():
    """Making sure all legal chains are created given a starting Residue() mapper."""
    ALL_NAMES = list(string.ascii_uppercase)
    result1 = sp.legal_chain_names(dict())
    assert set(result1) == set(ALL_NAMES)

    dummy_mapper = dict(zip(ALL_NAMES, range(len(ALL_NAMES))))
    result2 = sp.legal_chain_names(dummy_mapper)
    assert not result2

    ALL_NAMES.remove("A")
    dummy_mapper = dict(zip(ALL_NAMES, range(len(ALL_NAMES))))
    result3 = sp.legal_chain_names(dummy_mapper)
    assert result3 == ["A"]


def test_name_chains():
    """Ensuring that the name_chains() correctly names new chains."""

    def get_chains(fname) -> Dict[str, Chain]:
        """Helper testing method to get the chains from a PDB file."""
        reader = PandasPdb()
        reader.read_pdb(fname)
        res_mapper: Dict[str, Residue] = sp.build_residues(reader.df["ATOM"])
        chain_mapper: Dict[str, Chain] = sp.build_chains(res_mapper)
        return chain_mapper

    two_chain = f"{DATA_DIR}/two_chain.pdb"
    three_chain = f"{DATA_DIR}/three_chain.pdb"
    four_chain = f"{DATA_DIR}/four_chain.pdb"
    two_mapper: Dict[str, Chain] = get_chains(two_chain)
    three_mapper: Dict[str, Chain] = get_chains(three_chain)
    four_mapper: Dict[str, Chain] = get_chains(four_chain)
    two_mapper = sp.name_chains(two_mapper)
    three_mapper = sp.name_chains(three_mapper)
    four_mapper = sp.name_chains(four_mapper)

    assert set(two_mapper.keys()) == {"A", "B"}
    assert len(two_mapper["A"]) == 2
    assert len(two_mapper["B"]) == 5
    assert set(three_mapper.keys()) == {"A", "B", "C"}
    assert len(two_mapper["A"]) == 2
    assert len(three_mapper["B"]) == 3
    assert len(three_mapper["C"]) == 2
    assert set(four_mapper.keys()) == {"A", "B", "C", "D"}
    assert len(four_mapper["A"]) == 2
    assert len(four_mapper["B"]) == 3
    assert len(four_mapper["C"]) == 1
    assert len(four_mapper["D"]) == 1


def test_name_chains_already_named():
    """Ensuring that the name_chains() does not changed already named chains."""
    already_named = {"A": [], "B": []}
    already_named = sp.name_chains(already_named)
    assert set(already_named.keys()) == {"A", "B"}


def test_categorize_residue_all_canonical():
    """Checking that the structure_parser.categorize_residue() works for only canonical Residues()."""
    pdb_file = f"{DATA_DIR}/two_chain.pdb"
    structure: Structure = structure_from_pdb(pdb_file)
    assert structure.num_chains == 2
    # TODO(CJ): maybe make this a function for the Structure() class
    residues: List[Residue] = structure.residues
    for rr in residues:
        new_rr = sp.categorize_residue(rr)
        assert isinstance(new_rr, Residue)
        assert new_rr.is_canonical()
        assert not new_rr.is_metal()
        assert not new_rr.is_ligand()
        assert new_rr.rtype() == chem.ResidueType.CANONICAL


def test_categorize_residue_metal():
    """Checking that the structure_parser.categorize_residue() works for something that should become a MetalAtom()."""
    warnings.filterwarnings("ignore")
    pdb_file = f"{DATA_DIR}/just_metal.pdb"
    reader = PandasPdb()
    reader.read_pdb(pdb_file)
    zn_atom = Atom(**reader.df["HETATM"].iloc[0])
    base_residue = Residue("A.ZN.500", [zn_atom])
    metal: MetalAtom = sp.categorize_residue(base_residue)
    assert isinstance(metal, MetalAtom)
    assert not metal.is_canonical()
    assert metal.is_metal()
    assert not metal.is_ligand()
    assert metal.rtype() == chem.ResidueType.METAL


def test_categorize_residue_ligand():
    """Checking that the structure_parser.categorize_residue() works for something that should become a Ligand()."""

    pdb_file = f"{DATA_DIR}/just_ligand.pdb"
    reader = PandasPdb()
    reader.read_pdb(pdb_file)
    atoms: List[Residue] = list(
        map(lambda pr: Atom(**pr[-1]), reader.df["ATOM"].iterrows())
    )
    base_residue = Residue(".FAH.1", atoms)
    ligand: Ligand = sp.categorize_residue(base_residue)
    assert isinstance(ligand, Ligand)
    assert not ligand.is_canonical()
    assert not ligand.is_metal()
    assert ligand.is_ligand()
    assert ligand.rtype() == chem.ResidueType.LIGAND


def test_categorize_residue_solvent():
    """Checking that the structure_parser.categorize_residue() works for something that should become a Solvent()."""
    base_residue = Residue("A.HOH.1", [Atom(line_idx=1), Atom(line_idx=2)])
    solvent: Solvent = sp.categorize_residue(base_residue)
    assert isinstance(solvent, Solvent)
    assert id(base_residue) != id(solvent)
    for a1, a2 in zip(base_residue.atoms, solvent.atoms):
        assert id(a1) != id(a2)
    assert solvent.is_rd_solvent()
    assert solvent.rtype() == chem.ResidueType.SOLVENT


def test_build_chains():
    """Ensuring the structure_parser.build_chains() method correctly aggregates Residue() objects into Chain() objects."""
    pdb_file = f"{DATA_DIR}/four_chain.pdb"
    reader = PandasPdb()
    reader.read_pdb(pdb_file)
    res_mapper: Dict[str, Residue] = sp.build_residues(
        pd.concat((reader.df["ATOM"], reader.df["HETATM"]))
    )
    chain_mapper: Dict[str, Chain] = sp.build_chains(res_mapper)

    assert set(chain_mapper.keys()) == {"A", "B", "C", "D"}
    assert len(chain_mapper["A"]) == 2
    print(chain_mapper["B"].residues())
    assert len(chain_mapper["B"]) == 3
    assert len(chain_mapper["C"]) == 1
    assert len(chain_mapper["D"]) == 1

    atom_counts = [12, 15, 7, 6]
    for ac, chain in zip(atom_counts, chain_mapper.values()):
        assert chain.num_atoms() == ac
        for res in chain.residues():
            assert isinstance(res, Residue)


def test_structure_from_pdb_bad_input():
    """Checking that the main method structure_parser.structure_from_pdb() fails for bad input."""
    with pytest.raises(SystemExit) as exe:
        sp.structure_from_pdb("dne.pdb")

    assert exe
    assert exe.type == SystemExit
    assert exe.value.code == 1


def test_structure_from_pdb_simple():
    """Checking that the main method structure_parser.structure_from_pdb() works for simple cases."""
    pdb_name = f"{DATA_DIR}/two_chain.pdb"
    structure: Structure = structure_from_pdb(pdb_name)

    chain_A_target = ["A.ALA.10", "A.THR.11"]
    chain_B_target = ["B.GLY.12", "B.GLY.13", "B.ASN.14", "B.LEU.15", "B.PRO.16"]
    assert structure.num_chains == 2
    (chain_A, chain_B) = structure.chains
    assert chain_A_target == list(map(str, chain_A))
    assert chain_B_target == list(map(str, chain_B))


def test_structure_from_pdb_mixed_case():
    """Checking that the main method structure_parser.structure_from_pdb() works for cases with mixed Residue() types.."""
    mixed_pdb = f"{DATA_DIR}/mixed_residues.pdb"
    structure: Structure = structure_from_pdb(mixed_pdb)
    assert structure.num_chains == 5
    all_residues = structure.residues

    assert all_residues[0].is_canonical()
    assert all_residues[1].is_canonical()
    assert all_residues[2].is_canonical()
    assert all_residues[3].is_canonical()
    assert all_residues[4].is_ligand()
    assert all_residues[5].is_metal()
    assert all_residues[6].is_rd_solvent()

    assert isinstance(all_residues[0], Residue)
    assert isinstance(all_residues[1], Residue)
    assert isinstance(all_residues[2], Residue)
    assert isinstance(all_residues[3], Residue)
    assert isinstance(all_residues[4], Ligand)
    assert isinstance(all_residues[5], MetalAtom)
    assert isinstance(all_residues[6], Solvent)


def test_ligand_from_pdb():
    """Checking that the method structure_parser.ligand_from_pdb() works for a simple case and returns a Ligand()."""
    ligand_pdb = f"{DATA_DIR}/just_ligand.pdb"
    ligand = sp.ligand_from_pdb(ligand_pdb)
    assert ligand.is_ligand()
    assert len(ligand.atoms) == 7
    assert ligand.name == "FAH"


def test_ligand_from_pdb_bad_input():
    """Checking that the method structure_parser.ligand_from_pdb() fails for a bad input."""
    with pytest.raises(SystemExit) as exe:
        sp.ligand_from_pdb("dne.pdb")

    assert exe
    assert exe.type == SystemExit
    assert exe.value.code == 1

#TODO(CJ): add tests for structure_parser.get_ligand_name()
