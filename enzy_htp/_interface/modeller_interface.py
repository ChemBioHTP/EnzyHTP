"""Defines a ModellerInterface class that serves as a bride for enzy_htp to utilize the Modeller software package.
Uses the ModellerConfig class found in enzy_htp/_config/modeller_config.py Supported operations include:

    + filling missing loop Residue()'s in a Structure()

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2024-05-30
"""
import os
from typing import List

from pathlib import Path

from .base_interface import BaseInterface

import importlib


from enzy_htp import config
from enzy_htp.chemical import SeqRes

from enzy_htp.core import file_system as fs
from enzy_htp.core.general import HiddenPrints

from enzy_htp._config.modeller_config import ModellerConfig, default_modeller_config

from enzy_htp.structure import (
    Chain,
    Structure,
    Residue,
    PDBParser
    )


class ModellerInterface(BaseInterface):
    """Class that provides a direct interface for enzy_htp to utilize Modeller. Supported operations
    include adding missing Residue()'s Users should use this class as the only way to interact with this
    Application.

    Attributes:
        config_ : The ModellerConfig() class which provides settings for both running Modeller and maintaining a compatible environment.
        env_manager_ : The EnvironmentManager() class which ensures all required environment elements exist.
        compatible_env_ : A bool() indicating if the current environment is comptaible with the object itself.
        self.modeller_ : The modeller package import that is accesible only within the class.
        self.modeller_automodel_ : The modeller.automodel package import that is accesible only within the class.
    """

    def __init__(self, parent, config: ModellerConfig = None) -> None:
        """Simplistic constructor that optionally takes an AlphaFillConfig object as its only argument.
        Calls parent constructor.
        """
        super().__init__(parent, config, default_modeller_config)
        self.modeller_ = None
        self.modeller_automodel_ = None
        try:
            self.modeller_ = importlib.import_module('modeller')
        except:
            pass

        try:
            self.modeller_automodel_ = importlib.import_module('modeller.automodel')
        except:
            pass


    def delete_temp_files(self, stem:str) -> None:
        """Removes various temporary files created by modelller during loop modelling."""
        for ext in """.BL00010001.pdb
        .B99990001.pdb
        .BL00020001.pdb
        .D00000001
        .DL00010001
        .DL00020001
        .IL00000001.pdb
        .ini
        .lrsr
        .rsr
        .sch
        .V99990001""".split():
            fs.safe_rm( f"{stem}{ext}" )


    @property
    def modeller(self) -> "module":
        """Gets the modeller module if it exists, raises an error if not."""
        if self.modeller_ is None:
            err_msg:str="The 'modeller' python package is not installed in this environment. Cannot use ModellerInterface()."
            _LOGGER.error( err_msg )  
            raise ImportError(err_msg)
        return self.modeller_


    @property
    def modeller_automodel(self) -> "module":
        """Gets the modelle._automodel module if it exists, raises an error if not."""
        if self.modeller_automodel_ is None:
            err_msg:str="The 'modeller_automodel' python package is not installed in this environment. Cannot use ModellerInterface()."
            _LOGGER.error( err_msg )  
            raise ImportError(err_msg)
        return self.modeller_automodel_


    def add_missing_residues(self,
                            stru:Structure,
                            missing_residues:List[SeqRes],
                            work_dir:str=None,
                            **kwargs) -> None:
        """Uses the LoopMode class to add missing loop residues to a Structure(). All work is done 
        to the Structure() in place. Note that it will implicitly relax most/all of the Structure() including
        various non-loop sidechains.


        Args:
            stru: The Structure() to add missing Residue()'s to.
            missing_residues: The List[SeqRes] of missing residues to add.
            work_dir: The name of the directory where the work should be done. Optional.
        
        Returns:
            Nothing.
        """

        if work_dir is None:
            work_dir = config['system.SCRATCH_DIR']

        fs.safe_mkdir( work_dir )

        existing_residues = list()

        for rr in stru.residues:
            existing_residues.append( rr.create_seq_res() )


        all_residues = sorted(existing_residues + missing_residues)

        current_seq = str()
        target_seq = str()

        first_existing = None
        last_existing = None
    
        sidx:int=0
        for aidx, ar in enumerate(all_residues):
            
            if not ar.is_canonical():
                continue

            chain:Chain=stru.get_chain(ar.chain)
            if chain is None:
                continue
            
            if aidx > 0 and (ar.chain != all_residues[aidx-1].chain):
                current_seq += '/'                                
                target_seq += '/'                                
            
            ar.seq_idx = sidx
            sidx += 1
            target_seq += ar.one_letter()
            if ar.missing: 
                current_seq += '-' 
            else:
                if not first_existing:
                    first_existing = ar
                last_existing = ar
                current_seq += ar.one_letter()
        
        sp = PDBParser()

        current_seq += '*'
        target_seq += '*'

        align_file:str = f"{work_dir}/alignment.ali"
        temp_file:str = f"{work_dir}/modeller_temp.pdb"

        sp.save_structure( temp_file, stru )

        lines:List[str] = [
            ">P1;modeller_temp",
           f"structureX:modeller_temp:{first_existing.idx}:{first_existing.chain}:{last_existing.idx}:{last_existing.chain}:undefined:undefined:-1.00:-1.00",
           f"{current_seq}",
            ">P1;modeller_fill",
            "sequence:::::::::",
           f"{target_seq}"
        ]

        fs.write_lines(align_file, lines)
        #TODO(CJ): probably need to check for loops that are too long
        Environ = self.modeller.Environ
        LoopModel = self.modeller_automodel.LoopModel
        refine = self.modeller_automodel.refine
        log = self.modeller.log 

        log.none()
        env = Environ()

        start_dir:str=f"{Path(os.getcwd()).absolute()}"

        os.chdir(work_dir)
        
        # directories for input atom files
        env.io.atom_files_directory = ['.', '../atom_files']
        
        a = LoopModel(env, alnfile = 'alignment.ali',
                      knowns = 'modeller_temp', sequence = 'modeller_fill')
        a.starting_model= 1
        a.ending_model  = 1
        
        a.loop.starting_model = 1
        a.loop.ending_model   = 2
        a.loop.md_level       = refine.fast

        with HiddenPrints() as hp:
            a.make()
       
        fs.safe_mv("modeller_fill.B99990001.pdb", "modeller_fill.pdb")


        self.delete_temp_files( "modeller_fill" )

        os.chdir(start_dir)

        session = self.parent().pymol.new_session() 
        self.parent().pymol.general_cmd(session, [
            ('load', f"{work_dir}/modeller_fill.pdb"),
            ("load", f"{work_dir}/modeller_temp.pdb"),
            ('align', 'modeller_fill', 'modeller_temp'),
            ('save', f"{work_dir}/modeller_fill.pdb", 'modeller_fill')
        ])

        filled_stru = sp.get_structure(f"{work_dir}/modeller_fill.pdb")

        for ar in all_residues:
            if ar.seq_idx is None:
                continue
            fs_res:Residue = filled_stru.residues[ar.seq_idx]
            fs_res.idx = ar.idx

            key:str=f"{ar.chain}.{ar.idx}"

            if not stru.has_residue(key):
                chain:Chain=stru.get_chain(ar.chain)
                if chain is not None:
                    chain.add(fs_res, sort=False)

        for chain in stru.chains:
            chain.sort_residues()

