"""Defines an AlphafoldInterface class that serves as a bridge for enzy_htp to utilize Alphafold software.

Author: Gemini
Date: 2025-08-16
"""
from __future__ import annotations
import tempfile
from typing import Union, List, Optional, Dict
from pathlib import Path

from .base_interface import BaseInterface
from enzy_htp.structure import Structure
from enzy_htp.structure.structure_io.pdb_io import PDBParser
from enzy_htp._config.alphafold_config import AlphafoldConfig
from enzy_htp.core.job_manager import ClusterJob, ClusterJobConfig

class AlphafoldInterface(BaseInterface):
    """Class that provides a direct inteface for enzy_htp to utilize Alphafold software.
    """

    def __init__(self, parent, config: AlphafoldConfig = None) -> None:
        """Simplistic constructor that optionally takes an AlphafoldConfig object as its only argument.
        Calls parent class."""
        super().__init__(parent, config, AlphafoldConfig)

    def predict(self, sequences: List[str], cluster_job_config: Optional[Union[ClusterJobConfig, Dict]] = None, **kwargs) -> Structure:
        """Wrapper for running AlphaFold.

        This function is a wrapper for the `run` method. It handles the submission of the job to a cluster if
        `cluster_job_config` is provided.
        """
        # TODO make fasta file first

        if cluster_job_config:
            # Run on cluster (run as an array and use wait_to_array_end_plus)
            # TODO
            pass
        else:
            # Run locally
            pass

    def run(
        self,
        fasta_path: str,
        out_dir: Union[str, Path],
        num_models: int,
        num_recycles: int,
        database_source: Optional[str] = None,
        num_relax: int = 0,
        relax_max_iteration: int = 200,
        use_templates: bool = False,
        max_template_date: Optional[str] = None,
        model_preset: Optional[str] = "monomer_ptm",
        db_preset: str = "reduced_dbs",
        addition_options: Optional[List[str]] = None,
    ) -> Dict[str, str]:
        """TODO complete docstring"""
        out_dir = Path(out_dir)
        out_dir.mkdir(exist_ok=True)

        # TODO
        pass

    def make_job(
        self,
        fasta_path: str,
        out_dir: Union[str, Path],
        num_models: int,
        num_recycles: int,
        database_source: Optional[str] = None,
        num_relax: int = 0,
        relax_max_iteration: int = 200,
        use_templates: bool = False,
        max_template_date: Optional[str] = None,
        model_preset: Optional[str] = "monomer_ptm",
        db_preset: str = "reduced_dbs",
        addition_options: Optional[List[str]] = None,
    ) -> ClusterJob:
        """Create a cluster job for running AlphaFold."""
        # TODO