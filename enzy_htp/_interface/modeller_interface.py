#TODO(CJ): documentation
import os
from typing import List

from pathlib import Path

from .base_interface import BaseInterface

from enzy_htp import config

from enzy_htp.core import file_system as fs

from enzy_htp._config.modeller_config import ModellerConfig, default_modeller_config

from modeller import *
from modeller.automodel import *

from enzy_htp.structure import (
    Structure,
    Residue,
    PDBParser
    )


class ModellerInterface(BaseInterface):


    def __init__(self, parent, config: ModellerConfig = None) -> None:
        """Simplistic constructor that optionally takes an AlphaFillConfig object as its only argument.
        Calls parent constructor.
        """
        super().__init__(parent, config, default_modeller_config)



    def add_missing_residues(self,  stru, missing_residues, work_dir:str=None, inplace:bool=True):

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
        
        a.make()
       
        fs.safe_mv("modeller_fill.B99990001.pdb", "modeller_fill.pdb")

        for tk in """modeller_fill.BL00010001.pdb
        modeller_fill.B99990001.pdb
        modeller_fill.BL00020001.pdb
        modeller_fill.D00000001
        modeller_fill.DL00010001
        modeller_fill.DL00020001
        modeller_fill.IL00000001.pdb
        modeller_fill.ini
        modeller_fill.lrsr
        modeller_fill.rsr
        modeller_fill.sch
        modeller_fill.V99990001""".split():
            fs.safe_rm(tk)
   
        os.chdir(start_dir)

        session = self.parent().pymol.new_session() 
        self.parent().pymol.general_cmd(session, [
            ('load', f"{work_dir}/modeller_fill.pdb"),
            ("load", f"{work_dir}/modeller_temp.pdb"),
            ('align', 'modeller_fill', 'modeller_temp'),
            ('save', f"{work_dir}/modeller_fill.pdb", 'modeller_fill')
        ])

        filled_stru = sp.get_structure(f"{work_dir}/modeller_fill.pdb")

        if not inplace:
            temp = stru.clone()
            stru = temp

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


        if not inplace:
            return stru
