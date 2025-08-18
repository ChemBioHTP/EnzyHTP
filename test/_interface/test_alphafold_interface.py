"""Test module for enzy_htp._interface.alphafold_interface

Author: Gemini
Date: 2025-08-16
"""
import pytest
from unittest.mock import MagicMock, patch, ANY
from pathlib import Path

from enzy_htp._interface.alphafold_interface import AlphafoldInterface, AlphafoldConfig
from enzy_htp.core.job_manager import ClusterJobConfig
from enzy_htp.structure.structure_io.pdb_io import PDBParser
from enzy_htp.structure.structure import Structure


# TODO make unit test
# use "MSTPSLIPSGVHEVLAKYKDGN" for all the tests.