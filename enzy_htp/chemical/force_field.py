"""This module contains force field information.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2025-08-06
"""

AMBER_PROTEIN_FF_BACKBONE_ATOM_TYPE_MAPPER = {
    'fb15':         {'C': 'C','CA': 'CT','H': 'H','N': 'N','O': 'O' ,'OXT': 'O2'},
    'ff03.r1':      {'C': 'C','CA': 'CT','H': 'H','N': 'N','O': 'O' ,'OXT': 'O2'},
    'ff03ua':       {'C': 'C','CA': 'CT','H': 'H','N': 'N','O': 'O' ,'OXT': 'O2'},
    'ff14SB':       {'C': 'C','CA': 'CX','H': 'H','N': 'N','O': 'O' ,'OXT': 'O2'},
    'ff14SBonlysc': {'C': 'C','CA': 'CX','H': 'H','N': 'N','O': 'O' ,'OXT': 'O2'},
    'ff15ipq':      {'C': 'C','CA': 'CX','H': 'H','N': 'N','O': 'OD','OXT': 'O3'},
    'ff15ipq-vac':  {'C': 'C','CA': 'CX','H': 'H','N': 'N','O': 'OD','OXT': 'O3'},
    'ff19SB':       {'C': 'C','CA': 'XC','H': 'H','N': 'N','O': 'O' ,'OXT': 'O2'},
    'ff19ipq':      {'C': 'C','CA': 'CX','H': 'H','N': 'N','O': 'OD','OXT': 'O3'},
}
"""Dictionary mapping Amber22 force field names to backbone atom type mappings.
Generated from Amber22's leaprc files, this mapping provides for each force
field a nested dictionary that maps backbone atom names ('C', 'CA', 'H', 'N', 'O', 'OXT')
to their corresponding Amber atom types.
"""