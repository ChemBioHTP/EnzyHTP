"""Testing the enzy_htp.mutation.mutation_pattern.position_pattern.py submodule.
Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-02-17
"""

import logging
import os

from enzy_htp import PDBParser
from enzy_htp.core.logger import _LOGGER
import enzy_htp.mutation.mutation_pattern.position_pattern as m_p

CURR_FILE = os.path.abspath(__file__)
CURR_DIR = os.path.dirname(CURR_FILE)
DATA_DIR = f"{CURR_DIR}/../data/"
sp = PDBParser()


def test_decode_position_pattern():
    """test the function use a made up position_pattern for KE"""
    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_pattern = "resi 254 around 4"

    assert set(m_p.decode_position_pattern(test_stru, test_pattern)) == set([('A', 224), ('A', 202), ('A', 9), ('A', 201), ('A', 48),
                                                                             ('A', 11), ('A', 169), ('A', 222), ('A', 101), ('A', 128),
                                                                             ('A', 50)])


def test_decode_position_pattern_dimer():
    """test the function use a made up position_pattern for puo"""
    test_pdb = f"{DATA_DIR}puo_put.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_pattern = "resi 267 around 5"
    assert set(m_p.decode_position_pattern(test_stru, test_pattern)) == set([('A', 273), ('A', 268), ('A', 263), ('A', 389), ('A', 290),
                                                                             ('A', 269), ('A', 264), ('A', 390), ('A', 270), ('A', 265),
                                                                             ('A', 391), ('A', 271), ('A', 266), ('A', 388),
                                                                             ('B', 719)])  # note ligand wont be included

def test_decode_builtin_function():
    """test the function use a made up builtin_function pattern for KE"""

    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_pattern = "$ef_hotspot('B.254.CAE', 'B.254.H2', (170,180))"

    existing_level = _LOGGER.level
    _LOGGER.setLevel(logging.DEBUG)

    result = m_p.decode_builtin_function(test_stru, test_pattern)

    _LOGGER.setLevel(existing_level)

    assert set([res.key_str for res in result]) == set(["A.21", "A.24", "A.22", "A.228"])

def test_decode_builtin_function_vec():
    """test the function use a made up builtin_function pattern for KE.
    test the dispatch on using vector"""
    existing_level = _LOGGER.level
    _LOGGER.setLevel(logging.DEBUG)

    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_pattern = "$ef_hotspot((-0.373, -0.1285, 0.369), (22.665,-2.5785,-52.461), (170,180))"
    result = m_p.decode_builtin_function(test_stru, test_pattern)
    _LOGGER.setLevel(existing_level)
    assert set([res.key_str for res in result]) == set(["A.21", "A.24", "A.22", "A.228"])

def test_decode_builtin_function_surface():
    """test the function use a made up builtin_function pattern for KE."""

    test_pdb = f"{DATA_DIR}KE_07_R7_2_S.pdb"
    test_stru = sp.get_structure(test_pdb)
    test_pattern = "$surface()"
    result = m_p.decode_builtin_function(test_stru, test_pattern)
    assert set([res.key_str for res in result]) == {'A.1', 'A.104', 'A.105', 'A.108', 'A.109', 'A.11', 'A.110', 'A.111', 'A.114', 'A.115', 'A.118', 'A.119', 'A.120', 'A.121', 'A.122', 'A.123', 'A.13', 'A.132', 'A.133', 'A.134', 'A.135', 'A.136', 'A.137', 'A.14', 'A.143', 'A.144', 'A.145', 'A.146', 'A.147', 'A.148', 'A.149', 'A.15', 'A.150', 'A.151', 'A.152', 'A.154', 'A.155', 'A.158', 'A.159', 'A.16', 'A.161', 'A.162', 'A.163', 'A.164', 'A.166', 'A.173', 'A.174', 'A.175', 'A.176', 'A.177', 'A.178', 'A.179', 'A.18', 'A.180', 'A.184', 'A.185', 'A.188', 'A.189', 'A.19', 'A.191', 'A.192', 'A.193', 'A.195', 'A.196', 'A.197', 'A.2', 'A.20', 'A.202', 'A.205', 'A.206', 'A.207', 'A.208', 'A.21', 'A.211', 'A.212', 'A.215', 'A.217', 'A.219', 'A.22', 'A.224', 'A.225', 'A.229', 'A.230', 'A.231', 'A.233', 'A.235', 'A.236', 'A.239', 'A.24', 'A.240', 'A.242', 'A.243', 'A.244', 'A.245', 'A.247', 'A.249', 'A.25', 'A.251', 'A.252', 'A.253', 'A.26', 'A.27', 'A.28', 'A.29', 'A.3', 'A.31', 'A.34', 'A.37', 'A.38', 'A.40', 'A.41', 'A.42', 'A.45', 'A.46', 'A.5', 'A.52', 'A.54', 'A.55', 'A.56', 'A.57', 'A.58', 'A.59', 'A.60', 'A.61', 'A.64', 'A.67', 'A.68', 'A.70', 'A.71', 'A.72', 'A.74', 'A.75', 'A.76', 'A.82', 'A.84', 'A.85', 'A.86', 'A.87', 'A.91', 'A.93', 'A.94', 'A.95', 'A.96', 'A.98'}
