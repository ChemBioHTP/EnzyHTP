"""Testing the PrmtopParser class in the enzy_htp.structure.structure_io.prmtop_io

Author: QZ Shao <shaoqz@icloud.com>
Date: 2024-01-23
"""
from pathlib import Path

import pytest

import enzy_htp.core.file_system as fs
from enzy_htp.structure.structure_io import PrmtopParser

BASE_DIR = Path(__file__).absolute().parent
DATA_DIR = f"{BASE_DIR}/../data"

def test_prmtop_parser_get_stru():
    """make sure function works as expected"""
    test_prmtop = f"{DATA_DIR}/KE_07_R7_2_S_sol_addpdb.prmtop"
    test_stru = PrmtopParser().get_structure(test_prmtop)
    assert test_stru.residues[0].num_atoms == 19
    assert test_stru.residues[0].find_atom_name("N").charge == 0.1592
    assert test_stru.residues[-1].idx == 6700
    assert test_stru.residues[-1].name == "WAT"
    assert set(test_stru.residues[253].atom_name_list) == set([
        "OAB","NAL","OAA","CAI","CAF","H1",
        "CAC","H4","CAD","H3","CAK","CAJ",
        "CAE","H2","NAG","OAH",
    ])
    assert test_stru.atoms[-1].charge == 0.417

def test_parse_version():
    """test using content from an example file: {DATA_DIR}/KE_07_R7_2_S_10f.prmtop"""
    test_content = "%VERSION  VERSION_STAMP = V0001.000  DATE = 01/03/24 16:27:03                  \n"
    assert PrmtopParser._parse_version(test_content) == {
        "VERSION_STAMP" : "V0001.000",
        "DATE" : "01/03/24 16:27:03",
    }

def test_parse_title():
    """test using content from an example file: {DATA_DIR}/KE_07_R7_2_S_10f.prmtop"""
    test_content = " TITLE                                                                     \n%FORMAT(20a4)                                                                   \ndefault_name                                                                    \n"
    assert PrmtopParser._parse_title(test_content) == {
        "TITLE" : "default_name"
    }

def test_parse_atom_name():
    """test using content from an example file: {DATA_DIR}/KE_07_R7_2_S_10f.prmtop"""
    test_content = """ ATOM_NAME                                                                 
%FORMAT(20a4)                                                                   
N   H1  H2  H3  CA  HA  CB  HB2 HB3 CG  HG2 HG3 SD  CE  HE1 HE2 HE3 C   O   N   
H   CA  HA  CB  HB2 HB3 CG  HG  CD1 HD11HD12HD13CD2 HD21HD22HD23C   O   N   H   
CA  HA  CB  HB1 HB2 HB3 C   O   N   H   CA  HA  CB  HB2 HB3 CG  HG2 HG3 CD  HD2 
HD3 CE  HE2 HE3 NZ  HZ1 HZ2 HZ3 C   O   N   H   CA  HA  CB  HB2 HB3 CG  HG2 HG3 
CD  HD2 HD3 NE  HE  CZ  NH1 HH11HH12NH2 HH21HH22C   O   N   H   CA  HA  CB  HB  
CG2 HG21HG22HG23CG1 HG12HG13CD1 HD11HD12HD13C   O   N   H   CA  HA  CB  HB2 HB3 
CG  OD1 OD2 C   O   N   H   CA  HA  CB  HB1 HB2 HB3 C   O   N   H   CA  HA  CB  
HB1 HB2 HB3 C   O   N   H   CA  HA  CB  HB2 HB3 CG  HG  CD1 HD11HD12HD13CD2 HD21
HD22HD23C   O   N   H   CA  HA  CB  HB  CG2 HG21HG22HG23CG1 HG12HG13CD1 HD11HD12
HD13C   O   N   H   CA  HA  CB  HB2 HB3 CG  HG2 HG3 SD  CE  HE1 HE2 HE3 C   O   
N   H   CA  HA  CB  HB2 HB3 CG  HG2 HG3 CD  HD2 HD3 CE  HE2 HE3 NZ  HZ1 HZ2 HZ3 
C   O   N   H   CA  HA  CB  HB2 HB3 CG  OD1 OD2 C   O   N   H   CA  HA2         """

    assert len(PrmtopParser._parse_atom_name(test_content)["ATOM_NAME"]) == 238

def test_parse_residue_number():
    """test using content from an example file: {DATA_DIR}/KE_07_R7_2_S_10f.prmtop
    after add_pdb"""
    test_content = """ RESIDUE_NUMBER
%COMMENT Residue number (resSeq) read from PDB file; DIMENSION(NRES)
%COMMENT Residue insertion code (iCode) not present in PDB file
%COMMENT If present, look for: %FLAG RESIDUE_ICODE
%FORMAT(20I4)
   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40
  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57  58  59  60
  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79  80
  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99 100
 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120
 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140
 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160
 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180
 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200
 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220
 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240
 241 242 243 244 245 246 247 248 249 250 251 252 253 254"""

    assert len(PrmtopParser._parse_residue_number(test_content)["RESIDUE_NUMBER"]) == 254

def test_parse_prmtop_file():
    """test using an example file"""
    data = PrmtopParser._parse_prmtop_file(f"{DATA_DIR}/prmtop_1")

    assert data

    n_atoms: int = 23231
    num_bond: int = 85
    num_ang: int = 186
    n_ptra: int = 190
    n_types: int = 19
    n_bond_h: int = 21805
    n_bond_a: int = 1450
    n_thet: int = 3319
    n_phi: int = 6138
    n_phi_h: int = 6551
    nnb: int = 43036
    n_phb: int = 1
    n_res: int = 6967

    assert len(data['ATOM_NAME']) == n_atoms
    assert len(data['CHARGE']) == n_atoms
    assert len(data['ATOMIC_NUMBER']) == n_atoms
    # assert len(data['MASS']) == n_atoms NOTE add these when they are needed.
    # assert len(data['ATOM_TYPE_INDEX']) == n_atoms
    # assert len(data['NUMBER_EXCLUDED_ATOMS']) == n_atoms
    # assert len(data['NONBONDED_PARM_INDEX']) == n_types**2
    assert len(data['RESIDUE_LABEL']) == n_res
    assert len(data['RESIDUE_POINTER']) == n_res
    # assert len(data['BOND_FORCE_CONSTANT']) == num_bond
    # assert len(data['BOND_EQUIL_VALUE']) == num_bond
    # assert len(data['ANGLE_FORCE_CONSTANT']) == num_ang
    # assert len(data['ANGLE_EQUIL_VALUE']) == num_ang
    # assert len(data['DIHEDRAL_FORCE_CONSTANT']) == n_ptra
    # assert len(data['DIHEDRAL_PERIODICITY']) == n_ptra
    # assert len(data['DIHEDRAL_PHASE']) == n_ptra
    # assert len(data['SCEE_SCALE_FACTOR']) == n_ptra
    # assert len(data['SOLTY']) == 52
    # assert len(data['LENNARD_JONES_ACOEF']) == (n_types * (n_types + 1)) / 2
    # assert len(data['LENNARD_JONES_BCOEF']) == (n_types * (n_types + 1)) / 2
    # assert len(data['BONDS_INC_HYDROGEN']) == 3 * n_bond_h
    # assert len(data['BONDS_WITHOUT_HYDROGEN']) == 3 * n_bond_a
    # assert len(data['ANGLES_INC_HYDROGEN']) == 4 * n_thet
    # assert len(data['DIHEDRALS_INC_HYDROGEN']) == 5 * n_phi_h
    # assert len(data['DIHEDRALS_WITHOUT_HYDROGEN']) == 5 * n_phi
    # assert len(data['EXCLUDED_ATOMS_LIST']) == nnb
    # assert len(data['HBOND_ACOEF']) == n_phb
    # assert len(data['HBOND_BCOEF']) == n_phb
    # assert len(data['HBCUT']) == n_phb
    # assert len(data['AMBER_ATOM_TYPE']) == n_atoms
    # assert len(data['TREE_CHAIN_CLASSIFICATION']) == n_atoms
    # assert len(data['JOIN_ARRAY']) == n_atoms
    # assert len(data['IROTAT']) == n_atoms
    # assert len(data['RADII']) == n_atoms

def test_parse_prmtop_no_file():
    """Checking that PrmtopParser._parse_prmtop_file() throws an error if the listed file does not exist."""

    dne = Path('./dne')

    assert not dne.exists()

    with pytest.raises(ValueError) as exe:
        PrmtopParser._parse_prmtop_file(dne)

    assert exe

def test_parse_prmtop_file_add_pdb():
    """test using an example file after add pdb.
    (use a dry pdb for a solvated prmtop)"""
    data = PrmtopParser._parse_prmtop_file(f"{DATA_DIR}/KE_07_R7_2_S_sol_addpdb.prmtop")

    assert data

    n_atoms: int = 23302
    n_res: int = 6700

    assert len(data['ATOM_NAME']) == n_atoms
    assert len(data['CHARGE']) == n_atoms
    assert len(data['ATOMIC_NUMBER']) == n_atoms
    assert len(data['RESIDUE_LABEL']) == n_res
    assert len(data['RESIDUE_POINTER']) == n_res
    assert len(data['ATOM_NUMBER']) == n_atoms
    assert len(data['ATOM_ELEMENT']) == n_atoms
    assert len(data['RESIDUE_CHAINID']) == n_res
    assert len(data['RESIDUE_NUMBER']) == n_res
