"""Testing the CIFParser class in the enzy_htp.structure.structure_io.cif_io
This class is for parsing in and out from mmCIF to Structure

Author: Zhong, Yinjie <yinjie.zhong@outlook.com>
Date: 2026-10-07
"""
import os

import numpy as np
import pandas as pd
import pytest
from biopandas.mmcif import PandasMmcif

from enzy_htp.core import file_system as fs
import enzy_htp.chemical as chem
from enzy_htp.structure.structure_io.cif_io import CIFParser, get_cif_df
from enzy_htp.structure import (
    Structure,
    Residue,
    Chain,
    Atom,
    MetalUnit,
    Ligand,
    Solvent,
    ModifiedResidue,
)

# pylint: disable=protected-access, invalid-name
CURRDIR = os.path.dirname(os.path.abspath(__file__))
WORK_DIR = f'{CURRDIR}/../work_dir/'
# the test data. 3FCR is a selenomethionine-labeled protein holding everything a
# mmCIF parser should handle: a single auth chain id (A) shared by the polymer,
# the modified residues (MSE x6, LLP x1 written as HETATM records sharing the
# polymer entity id), the PLP ligand, the EDO trash ligands, waters, and alt loc
# records (A/B/C)
CIF_FILE_3FCR = os.path.normpath(os.path.join(CURRDIR, '..', 'data', '3FCR.cif'))
# 1NVG is a real RCSB entry (a zinc-finger protein) used for the PDB/CIF parser
# consistency verification: waters and ZN ions share the auth chain A of the
# polymer and the HIS 39 side chain has alt loc A/B records
CIF_FILE_1NVG = os.path.normpath(os.path.join(CURRDIR, '..', 'data', '1NVG.cif'))
sp = CIFParser()

# a mini cif with alt loc records. residue GLY 1 has CA with alt loc A/B and
# residue ALA 2 has CB with only alt loc A
ALT_LOC_CIF = f"""data_test
#
_entry.id
test
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . GLY A 1 1 ? 0.0 0.0 0.0 1.00 10.00 ? 1 GLY A N 1
ATOM 2 CA CA A GLY A 1 1 ? 1.0 0.0 0.0 1.00 10.00 ? 1 GLY A CA 1
ATOM 3 CA CA B GLY A 1 1 ? 2.0 0.0 0.0 1.00 10.00 ? 1 GLY A CA 1
ATOM 4 C C . GLY A 1 1 ? 3.0 0.0 0.0 1.00 10.00 ? 1 GLY A C 1
ATOM 5 N N . ALA A 1 2 ? 4.0 0.0 0.0 1.00 10.00 ? 2 ALA A N 1
ATOM 6 CB CB A ALA A 1 2 ? 5.0 0.0 0.0 1.00 10.00 ? 2 ALA A CB 1
ATOM 7 O O . ALA A 1 2 ? 6.0 0.0 0.0 1.00 10.00 ? 2 ALA A O 1
#
"""

# a mini cif with 2 models of the same GLY residue at different coordinates
MULTI_MODEL_CIF = f"""data_test
#
_entry.id
test
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . GLY A 1 1 ? 0.0 0.0 0.0 1.00 10.00 ? 1 GLY A N 1
ATOM 2 CA CA . GLY A 1 1 ? 1.0 0.0 0.0 1.00 10.00 ? 1 GLY A CA 1
ATOM 3 C C . GLY A 1 1 ? 2.0 0.0 0.0 1.00 10.00 ? 1 GLY A C 1
ATOM 4 N N . GLY A 1 1 ? 5.0 0.0 0.0 1.00 10.00 ? 1 GLY A N 2
ATOM 5 CA CA . GLY A 1 1 ? 6.0 0.0 0.0 1.00 10.00 ? 1 GLY A CA 2
ATOM 6 C C . GLY A 1 1 ? 7.0 0.0 0.0 1.00 10.00 ? 1 GLY A C 2
#
"""

# a mini cif mimicking a file converted from the PDB format: the modified residue
# MSE is written as a HETATM record but shares the polymer entity id (1) with a
# numeric label_seq_id, while the ligand PLP, the metal CA, and the solvent HOH
# share the same auth chain id yet have their own entity ids and a blank label_seq_id
MODIFIED_CIF = f"""data_test
#
_entry.id
test
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . GLY A 1 1 ? 0.0 0.0 0.0 1.00 10.00 ? 1 GLY A N 1
ATOM 2 CA CA . GLY A 1 1 ? 1.0 0.0 0.0 1.00 10.00 ? 1 GLY A CA 1
ATOM 3 C C . GLY A 1 1 ? 2.0 0.0 0.0 1.00 10.00 ? 1 GLY A C 1
HETATM 4 SE SE . MSE A 1 2 ? 3.0 0.0 0.0 1.00 10.00 ? 2 MSE A SE 1
HETATM 5 CE CE . MSE A 1 2 ? 3.5 0.0 0.0 1.00 10.00 ? 2 MSE A CE 1
HETATM 6 P1 P1 . PLP A 2 . ? 4.0 0.0 0.0 1.00 10.00 ? 3 PLP A P1 1
HETATM 7 O1 O1 . PLP A 2 . ? 4.5 0.0 0.0 1.00 10.00 ? 3 PLP A O1 1
HETATM 8 CA CA . CA A 5 . ? 5.5 0.0 0.0 1.00 10.00 ? 4 CA A CA 1
HETATM 9 O O . HOH A 4 . ? 5.0 0.0 0.0 1.00 10.00 ? 5 HOH A O 1
#
"""


def write_mini_cif(content: str, fname: str) -> str:
    """Helper method to write a mini cif into the work dir for testing."""
    os.makedirs(WORK_DIR, exist_ok=True)
    cif_path = f'{WORK_DIR}{fname}'
    fs.write_lines(cif_path, content.splitlines())
    return cif_path


def test_check_valid_cif_good_input():
    '''Good input for the check_valid_cif() helper method.'''
    # make dummy file
    dummy_cif = f'{CURRDIR}/dummy.cif'
    assert not os.path.exists(dummy_cif)
    fs.write_lines(dummy_cif, ['line1', 'line2'])
    # test
    assert not sp._check_valid_cif(dummy_cif)
    fs.safe_rm(dummy_cif)
    assert not os.path.exists(dummy_cif)


def test_check_valid_cif_bad_input():
    '''Bad inputs for the check_valid_cif() helper method.'''
    # non cif file
    txt_file = 'not_cif.txt'
    with pytest.raises(ValueError) as exe:
        sp._check_valid_cif(txt_file)
    assert exe
    assert exe.type == ValueError

    # cif that doesn't exist
    cif_imaginary = 'not_real.cif'
    with pytest.raises(ValueError) as exe:
        sp._check_valid_cif(cif_imaginary)
    assert exe
    assert exe.type == ValueError

    # cif with non-ASCII text
    non_ascii_cif = f'{CURRDIR}/bad_cif.cif'
    assert not os.path.exists(non_ascii_cif)
    fs.write_lines(non_ascii_cif, ['日本人 中國的'])
    assert os.path.exists(non_ascii_cif)

    with pytest.raises(ValueError) as exe:
        sp._check_valid_cif(non_ascii_cif)
    assert exe
    assert exe.type == ValueError

    fs.safe_rm(non_ascii_cif)
    assert not os.path.exists(non_ascii_cif)


def test_get_target_model():
    '''test if _get_target_model is getting the correct model'''
    cif_file = CIF_FILE_3FCR
    input_cif = PandasMmcif()
    input_cif.read_mmcif(cif_file)
    combined_df = pd.concat((input_cif.df['ATOM'], input_cif.df['HETATM']), ignore_index=True)

    target_model_df = sp._get_target_model(combined_df, model=0)
    assert len(target_model_df) == 3917
    with pytest.raises(ValueError) as exe:
        sp._get_target_model(combined_df, model=1)
    assert exe
    assert exe.type == ValueError

    # multiple model case
    mdl_cif = write_mini_cif(MULTI_MODEL_CIF, 'mini_two_model.cif')
    input_cif = PandasMmcif()
    input_cif.read_mmcif(mdl_cif)
    combined_df = pd.concat((input_cif.df['ATOM'], input_cif.df['HETATM']), ignore_index=True)
    mdl_0_df = sp._get_target_model(combined_df, model=0)
    mdl_1_df = sp._get_target_model(combined_df, model=1)
    assert list(mdl_0_df['Cartn_x']) == [0.0, 1.0, 2.0]
    assert list(mdl_1_df['Cartn_x']) == [5.0, 6.0, 7.0]
    fs.safe_rm(mdl_cif)


def test_get_cif_df():
    '''test the normalized df has the biopandas pdb naming and the CIF
    placeholder values are cleaned'''
    df = get_cif_df(CIF_FILE_3FCR)
    assert len(df) == 3917
    assert set(df.columns) >= {
        'record_name', 'atom_number', 'atom_name', 'residue_name', 'chain_id',
        'residue_number', 'x_coord', 'y_coord', 'z_coord', 'occupancy',
        'b_factor', 'alt_loc', 'element_symbol', 'charge', 'line_idx'
    }
    assert set(df['chain_id'].unique()) == {'A'}
    assert list(df['line_idx']) == list(range(len(df)))
    # CIF placeholder ('.', '?') and the 'None' string converted by biopandas
    # are all cleaned
    assert (df['alt_loc'] == 'None').sum() == 0
    assert (df['alt_loc'] == '.').sum() == 0


@pytest.mark.interface
def test_get_structure_bad_input():
    with pytest.raises(ValueError) as exe:
        sp.get_structure('dne.cif')
    assert exe
    assert exe.type == ValueError


@pytest.mark.interface
def test_get_structure_simple():
    '''Checking that the main method get_structure() works for a simple case.
    3FCR has a single auth chain A holding the polymer, the modified residues,
    the PLP ligand and the waters sharing the same chain id (mmCIF allowed).'''
    cif_file = CIF_FILE_3FCR
    stru: Structure = sp.get_structure(cif_file)

    assert stru.num_chains == 1
    assert set(c.name for c in stru.chains) == {'A'}
    # 3917 raw rows - 68 alt loc duplicate rows - 5 EDO trash residues (20 atoms)
    assert stru.num_atoms == 3829
    assert stru['A'][0].name == 'GLY'
    assert stru['A'][0].idx == 0
    assert stru['A'][0].is_canonical()
    assert stru.get('A.458').name == 'PLP'
    assert len(list(stru['A'])) == 772


@pytest.mark.interface
def test_get_structure_mixed_case():
    '''Checking that the residues in the peptide chain are correctly categorized
    even when non-polymer components share the same chain id (mmCIF allowed).
    3FCR also contains modified residues (MSE x6, LLP x1) written as HETATM
    records sharing the polymer entity id (the PDB-converted convention).'''
    cif_file = CIF_FILE_3FCR
    stru: Structure = sp.get_structure(cif_file)
    all_residues = stru.residues

    assert len(list(filter(lambda x: x.is_canonical(), all_residues))) == 451
    assert len(list(filter(lambda x: x.is_modified(), all_residues))) == 7
    assert len(list(filter(lambda x: x.is_ligand(), all_residues))) == 1
    assert len(list(filter(lambda x: x.is_solvent(), all_residues))) == 313

    chain_A = stru['A']
    assert isinstance(chain_A[0], Residue)
    assert isinstance(stru.get('A.1'), ModifiedResidue)
    assert stru.get('A.1').name == 'MSE'
    assert isinstance(stru.get('A.288'), ModifiedResidue)
    assert stru.get('A.288').name == 'LLP'
    assert isinstance(stru.get('A.458'), Ligand)
    assert isinstance(stru.get('A.464'), Solvent)


@pytest.mark.interface
def test_get_structure_keep_trash():
    '''the EDO trash ligand should only be removed when remove_trash=True'''
    cif_file = CIF_FILE_3FCR
    stru: Structure = sp.get_structure(cif_file, remove_trash=False)
    assert stru.num_atoms == 3849
    trash = list(filter(lambda x: x.rtype is chem.ResidueType.TRASH, stru.residues))
    assert len(trash) == 5
    assert all(t.name == 'EDO' for t in trash)
    assert trash[0].chain.name == 'A'


@pytest.mark.interface
def test_get_structure_give_idx_map():
    '''the idx change mapper is always empty for mmCIF since the chain id is
    carried by every atom natively'''
    cif_file = CIF_FILE_3FCR
    stru, idx_change_mapper = sp.get_structure(cif_file, give_idx_map=True)
    assert isinstance(stru, Structure)
    assert idx_change_mapper == {}


@pytest.mark.interface
def test_get_structure_model_selection():
    '''test get_structure() selects the correct model from a multi-model file'''
    mdl_cif = write_mini_cif(MULTI_MODEL_CIF, 'mini_two_model.cif')

    stru_0: Structure = sp.get_structure(mdl_cif, model=0)
    assert stru_0.num_atoms == 3
    assert list(stru_0.atoms[0].coord) == [0.0, 0.0, 0.0]
    assert list(stru_0.atoms[1].coord) == [1.0, 0.0, 0.0]

    stru_1: Structure = sp.get_structure(mdl_cif, model=1)
    assert stru_1.num_atoms == 3
    assert list(stru_1.atoms[0].coord) == [5.0, 0.0, 0.0]
    assert list(stru_1.atoms[1].coord) == [6.0, 0.0, 0.0]

    with pytest.raises(ValueError) as exe:
        sp.get_structure(mdl_cif, model=2)
    assert exe
    assert exe.type == ValueError

    fs.safe_rm(mdl_cif)


@pytest.mark.interface
def test_get_structure_alt_loc():
    '''test get_structure() on a mmCIF that contains alt loc.'''
    alt_cif = write_mini_cif(ALT_LOC_CIF, 'mini_alt_loc.cif')
    # default: keep the lexicographically smallest alt loc
    stru: Structure = sp.get_structure(alt_cif)
    res_1 = stru.get('A.1')
    ca_atom = list(filter(lambda x: x.name == 'CA', res_1.atoms))[0]
    assert len(res_1.atoms) == 3
    assert ca_atom.alt_loc == 'A'
    assert list(ca_atom.coord) == [1.0, 0.0, 0.0]
    # the single alt loc residue is cleaned by default
    res_2 = stru.get('A.2')
    assert 'CB' in [a.name for a in res_2.atoms]
    assert list(filter(lambda x: x.name == 'CB', res_2.atoms))[0].alt_loc is None

    # keep a specific alt loc. note the single alt loc residue is governed by
    # the single_alt_loc strategy and is cleaned by default
    stru: Structure = sp.get_structure(alt_cif, alt_loc_keep='B')
    res_1 = stru.get('A.1')
    ca_atom = list(filter(lambda x: x.name == 'CA', res_1.atoms))[0]
    assert ca_atom.alt_loc == 'B'
    assert list(ca_atom.coord) == [2.0, 0.0, 0.0]
    res_2 = stru.get('A.2')
    assert 'CB' in [a.name for a in res_2.atoms]

    # keep a specific alt loc without a fallback. with the "keep" strategy for
    # the single alt loc residue, all its alt loc atoms are deleted since the
    # preferred alt loc does not exist in it
    stru: Structure = sp.get_structure(alt_cif, alt_loc_keep='B', single_alt_loc='keep')
    assert list(stru.get('A.1.CA').coord) == [2.0, 0.0, 0.0]
    assert 'CB' not in [a.name for a in stru.get('A.2').atoms]

    # keep a specific alt loc with a fallback
    stru: Structure = sp.get_structure(alt_cif, alt_loc_keep='B else A', single_alt_loc='keep')
    assert list(stru.get('A.1.CA').coord) == [2.0, 0.0, 0.0]
    # the fallback alt loc is used when the preferred one does not exist
    cb_atom = list(filter(lambda x: x.name == 'CB', stru.get('A.2').atoms))[0]
    assert cb_atom.alt_loc == 'A'

    # delete the single alt loc residue atoms
    stru: Structure = sp.get_structure(alt_cif, single_alt_loc='delete')
    res_1 = stru.get('A.1')
    assert len(res_1.atoms) == 3
    res_2 = stru.get('A.2')
    assert 'CB' not in [a.name for a in res_2.atoms]

    fs.safe_rm(alt_cif)


@pytest.mark.interface
def test_get_file_str():
    '''test for the getting cif file from `Structure` instance and the content
    structure of the file'''
    cif_file = CIF_FILE_3FCR
    stru: Structure = sp.get_structure(cif_file)
    result_str = sp.get_file_str(stru)
    lines = result_str.splitlines()

    assert lines[0] == f'data_{CIFParser.DATA_BLOCK_NAME}'
    assert 'loop_' in lines
    assert lines[-1] == '#'
    tags = list(filter(lambda x: x.startswith('_atom_site.'), lines))
    assert tags == [f'_atom_site.{tag}' for tag in CIFParser.ATOM_SITE_TAGS]
    data_rows = list(filter(lambda x: x.startswith(('ATOM', 'HETATM')), lines))
    assert len(data_rows) == 3829
    for row in data_rows[:10]:
        assert len(row.split()) == len(CIFParser.ATOM_SITE_TAGS)
    # canonical residues are ATOM records and non-polymer ones are HETATM
    assert data_rows[0].startswith('ATOM')
    plp_rows = list(filter(lambda x: x.split()[5] == 'PLP', data_rows))  # label_comp_id field
    assert plp_rows and all(row.startswith('HETATM') for row in plp_rows)
    # modified residues written as HETATM in the file go back to ATOM records
    mse_rows = list(filter(lambda x: x.split()[5] == 'MSE', data_rows))
    assert mse_rows and all(row.startswith('ATOM') for row in mse_rows)


@pytest.mark.interface
def test_get_file_str_not_renumber_atom():
    '''test if_renumber=False keeps the atom index'''
    cif_file = CIF_FILE_3FCR
    stru: Structure = sp.get_structure(cif_file)
    stru.atoms[100].idx = 84242

    result_str = sp.get_file_str(stru, if_renumber=False)
    data_rows = list(filter(lambda x: x.startswith(('ATOM', 'HETATM')), result_str.splitlines()))
    assert len(list(filter(lambda x: x.split()[1] == '84242', data_rows))) == 1

    result_str = sp.get_file_str(stru, if_renumber=True)
    assert '84242' not in result_str


@pytest.mark.interface
def test_get_file_str_dispatch_residue_atom():
    '''test the dispatch of get_file_str() with `Residue` and/or `Atom` instances only'''
    cif_file = CIF_FILE_3FCR
    stru: Structure = sp.get_structure(cif_file)
    res = stru['A'][0]

    result_str = sp.get_file_str(res, if_renumber=False)
    data_rows = list(filter(lambda x: x.startswith(('ATOM', 'HETATM')), result_str.splitlines()))
    assert len(data_rows) == len(res.atoms)

    result_str = sp.get_file_str(res.atoms[0])
    data_rows = list(filter(lambda x: x.startswith(('ATOM', 'HETATM')), result_str.splitlines()))
    assert len(data_rows) == 1


@pytest.mark.interface
def test_save_structure():
    '''test save_structure() writes the same content as get_file_str()'''
    cif_file = CIF_FILE_3FCR
    result_file = f'{WORK_DIR}3FCR_write_result.cif'

    stru: Structure = sp.get_structure(cif_file)
    result_path: str = sp.save_structure(result_file, stru)
    assert result_path == result_file
    assert os.path.exists(result_file)
    assert fs.lines_from_file(result_file) == sp.get_file_str(stru).splitlines()
    fs.safe_rm(result_file)


@pytest.mark.interface
def test_write_read_round_trip():
    '''Ensuring that the `Structure` instance parsed from mmCIF can be written back out
    and parsed again in a round trip without losing information.'''
    cif_file = CIF_FILE_3FCR
    result_file = f'{WORK_DIR}3FCR_roundtrip_result.cif'

    stru_1: Structure = sp.get_structure(cif_file)
    sp.save_structure(result_file, stru_1)
    stru_2: Structure = sp.get_structure(result_file)

    assert set(c.name for c in stru_1.chains) == set(c.name for c in stru_2.chains)
    assert stru_1.num_atoms == stru_2.num_atoms
    res_mapper_1 = {(r.chain.name, r.idx, r.name): r for r in stru_1.residues}
    res_mapper_2 = {(r.chain.name, r.idx, r.name): r for r in stru_2.residues}
    assert set(res_mapper_1.keys()) == set(res_mapper_2.keys())
    for res_key, res_1 in res_mapper_1.items():
        atom_mapper_2 = {a.name: a for a in res_mapper_2[res_key].atoms}
        assert len(atom_mapper_2) == len(res_1.atoms)
        for atom in res_1.atoms:
            assert np.allclose(atom.coord, atom_mapper_2[atom.name].coord)
    fs.safe_rm(result_file)


def test_cif_quote():
    '''test the CIF string quoting of _cif_quote()'''
    assert CIFParser._cif_quote('CA') == 'CA'
    assert CIFParser._cif_quote("O5'") == '"O5\'"'
    assert CIFParser._cif_quote('a"b') == "'a\"b'"
    assert CIFParser._cif_quote('a b') == '"a b"'
    assert CIFParser._cif_quote('') == '""'


def test_get_structure_modified_residue_entity_signal():
    '''test that a HETATM residue sharing the polymer entity id with a numeric
    label_seq_id is categorized as a modified residue, while HETATM residues
    with their own entity id and a blank label_seq_id are not (mimics a file
    converted from the PDB format)'''
    cif_file = write_mini_cif(MODIFIED_CIF, 'modified_test.cif')
    stru: Structure = sp.get_structure(cif_file)

    assert stru.get('A.2.SE').residue.rtype == chem.ResidueType.MODIFIED
    assert isinstance(stru.get('A.2.SE').residue, ModifiedResidue)
    # the ligand, the metal, and the solvent share the auth chain id but have
    # their own entity ids and a blank label_seq_id
    assert stru.get('A.3.P1').residue.rtype == chem.ResidueType.LIGAND
    assert isinstance(stru.get('A.3.P1').residue, Ligand)
    assert stru.get('A.4.CA').residue.name == 'CA'
    assert isinstance(stru.get('A.4.CA').residue, MetalUnit)
    assert stru.get('A.5.O').residue.rtype == chem.ResidueType.SOLVENT
    assert isinstance(stru.get('A.5.O').residue, Solvent)
    # a file with placeholder entity ids carries no entity signal and falls back to
    # the record-name based categorization (biopandas requires the tags to exist)
    no_signal_lines = []
    for line in MODIFIED_CIF.splitlines():
        tokens = line.split()
        if tokens and tokens[0] in ('ATOM', 'HETATM'):
            tokens[7] = '.'  # label_entity_id
            line = ' '.join(tokens)
        no_signal_lines.append(line)
    no_signal_cif = write_mini_cif('\n'.join(no_signal_lines), 'modified_no_signal_test.cif')
    stru_no_signal: Structure = sp.get_structure(no_signal_cif)
    assert stru_no_signal.get('A.2.SE').residue.rtype == chem.ResidueType.LIGAND
    fs.safe_rm(cif_file)
    fs.safe_rm(no_signal_cif)


@pytest.mark.interface
def test_get_structure_1nvg_consistency():
    '''consistency of CIFParser against PDBParser verified on the 1NVG pair: all
    the 2616 atoms of the legacy de-watered PDB file are recovered from the mmCIF
    with identical coordinates and categorization (checked against PDBParser,
    see the CIFParser docstring), plus the mmCIF-only content (120 waters on the
    auth chain) and strictly more information (the kept alt loc designation and
    the real occupancy of the HIS 39 alternates that the PDB file stripped)'''
    stru: Structure = sp.get_structure(CIF_FILE_1NVG, remove_trash=False)

    # single chain A: the waters and the metals stay on the auth chain of the polymer
    assert stru.num_chains == 1
    assert set(c.name for c in stru.chains) == {'A'}
    # 2616 polymer atoms + 2 ZN + 120 waters
    assert stru.num_atoms == 2736
    # metals and solvent are correctly categorized inside the peptide chain
    assert isinstance(stru.get('A.400'), MetalUnit)
    assert stru.get('A.400').name == 'ZN'
    assert isinstance(stru.get('A.500'), MetalUnit)
    assert isinstance(stru.get('A.501'), Solvent)
    assert stru.get('A.501').name == 'HOH'
    # the alt loc residue keeps the kept alternates' designation and real occupancy
    his_39 = stru.get('A.39')
    alt_atoms = [a for a in his_39.atoms if a.alt_loc]
    assert len(alt_atoms) == 6
    assert all(a.alt_loc == 'A' for a in alt_atoms)
    assert stru.get('A.39.CG').occupancy == pytest.approx(0.6)


def test_get_file_str_write_alt_loc():
    '''test that get_file_str() writes the alt loc of each atom instead of a hardcoded "."'''
    cif_file = write_mini_cif(ALT_LOC_CIF, 'alt_loc_write_test.cif')
    stru: Structure = sp.get_structure(cif_file)  # keep='first' keeps the alt loc A of GLY CA

    result_str = sp.get_file_str(stru, if_renumber=False)
    data_rows = list(filter(lambda x: x.startswith(('ATOM', 'HETATM')), result_str.splitlines()))
    ca_row = list(filter(lambda x: x.split()[3] == 'CA', data_rows))  # label_atom_id field
    assert len(ca_row) == 1
    assert ca_row[0].split()[4] == 'A'  # label_alt_id field
    # atoms without alt loc are written as "."
    n_row = list(filter(lambda x: x.split()[19] == 'N', data_rows))  # auth_atom_id field
    assert len(n_row) == 2
    assert all(row.split()[4] == '.' for row in n_row)
    fs.safe_rm(cif_file)
