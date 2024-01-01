"""Defines an XTBInterface class that serves as a bridge for enzy_htp to utilize the xtb
semi-empirical quantum mechanics (QM) package. Uses the XTBConfig class found in enzy_htp/_config/xtb_config.py.
Supported operations include:
    + single point energy calculations
    + geometry optimization 


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-08-17
"""
from pathlib import Path
from typing import Tuple, List

from .base_interface import BaseInterface

from enzy_htp import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp._config.xtb_config import XTBConfig, default_xtb_config
from enzy_htp.structure import Structure, Atom, StructureConstraint, PDBParser
from enzy_htp.structure.structure_region import StructureRegion, create_region_from_residues
from enzy_htp.structure.structure_constraint import (
    CartesianFreeze,
    DistanceConstraint,
    AngleConstraint,
    DihedralConstraint,
    ResiduePairConstraint
)




from plum import dispatch

class XTBInterface(BaseInterface):
    """Class that provdes a direct interface for enzy_htp to utilize the xtb semi empirical method. Supported methods
    include single energy calculations and geometry optimization
    
    Attributes:
        config_ : The XTBXonfig() class which provides settings for both running xtb and maintaing a compatible environment.
        env_manager_ : The EnvironmentManager() class which ensures all required environment elements exist.
        compatible_env_ : a bool() indicating if the current environment is compatible with the object itself.
    """

    def __init__(self, parent, config: XTBConfig = None) -> None:
        """Simplicstic constructor that requires the parent interface as an argument and optionally takes an XTBConfig instance.
        Calls parent constructor.
        """
        super().__init__(parent, config, default_xtb_config)

    def single_point(self, stru: Structure,
                            charge: int = 0,
                            spin: int = 1,
                            n_iter: int = -1,
                            n_proc: int = -1,
                            sele_str:str=None,
                            work_dir:str=None,
                            cterm_cap:str=None,
                            nterm_cap:str=None,
                            **kwargs
                            ) -> float:
        """Performs a single point energy calculation on the supplied file. Checks if the file exists and is one of the correct
        file formats. Errors if not. Returns a value in Hartrees. 
        #TODO(CJ): update this stuff
        Args:
            fname: Name of the input file to use.
            charge: Optional charge of the system as an int(). Defaults to 0.
            spin: Optional spin of the system as an int(). Defaults to 0.
            n_iter: Number of SCF iterations. Optional, if not specified uses value in XTBConfig.
            n_proc: Number of processes to use in calculation. Optional, if not specified uses value in XTBConfig.
            work_dir:TODO(CJ)

        Returns:
            Single point energy value in Hartrees.           
        """
        if work_dir is None:
            work_dir = "./"

        if n_iter == -1:
            n_iter = self.config_.N_ITER

        if n_proc == -1:
            n_proc = self.config_.N_PROC

        fname:str=f"{work_dir}/__temp_xtb.xyz"
        lines:List[str] = ["", "Title goes here"]

        if sele_str is not None:
            session = self.parent().pymol.new_session()
            residue_list=self.parent().pymol.get_residue_list(session, stru, sele_str=sele_str, work_dir=work_dir)
            sr:StructureRegion = create_structure_region(stru, residue_list, nterm_cap=nterm_cap, cterm_cap=cterm_cap)
            
            lines[0] = str(len(sr.atoms))

            for aa in sr.atoms:
                (x,y,z) = aa.coord
                lines.append(f"{aa.element} {x} {y} {z}")
        else:
            lines[0] = str(len(stru.atoms))

            for aa in stru.atoms:
                (x,y,z) = aa.coord
                lines.append(f"{aa.element} {x} {y} {z}")

        fs.write_lines(fname, lines)

        results = self.env_manager_.run_command(
            self.config_.XTB_EXE,
            ["--chrg", str(charge), "--iterations",
             str(n_iter), "--parallel",
             str(n_proc), "--norestart", "--sp", fname])
        
        fs.safe_rm( fname )
        self._remove_temp_files(str(Path(fname).parent))

        for ll in results.stdout.splitlines():
            if ll.find('TOTAL ENERGY') != -1:
                return float(ll.split()[3])

        _LOGGER.error(f"ERROR: Single point energy calculation for file '{fname}' did not contain energy value. Exiting...")
        exit( 1 )

    def geo_opt(self,
                stru:Structure,
                charge: int = 0,
                spin: int = 1,
                n_iter: int = -1,
                n_proc: int = -1,
                constraints:List[StructureConstraint]=None,
                freeze_backbone:bool=True,
                freeze_hydrogens:bool=True,
                sele_str:str=None,
                work_dir:str=None,
                nterm_cap:str=None,
                cterm_cap:str=None,
                **kwargs
                ) -> float:
        """TODO(CJ)

        Args:

        Returns:
            
        """

        if n_iter == -1:
            n_iter = self.config_.N_ITER

        if n_proc == -1:
            n_proc = self.config_.N_PROC

        if work_dir is None:
            work_dir = "./"

        cmd_args:List[str] = [
                        "--chrg", str(charge),
                        "--iterations", str(n_iter),
                        "--parallel", str(n_proc),
                        "--norestart"
        ]
        input_lines:List[str]=list()

        _parser = PDBParser()
        fname:str = f"{work_dir}/__temp_xtb.pdb"

        freeze_indices:List[int] = list()
        constraint_lines:List[str] = list()

        if sele_str is not None:
            session = self.parent().pymol.new_session()
            residue_list=self.parent().pymol.get_residue_list(session, stru, sele_str=sele_str, work_dir=work_dir)
            sr:StructureRegion = create_structure_region(stru, residue_list, nterm_cap=nterm_cap, cterm_cap=cterm_cap)
            
            pdb_lines:List[str] = list()
            for aa in sr.atoms:
                pdb_lines.append( _parser._write_pdb_atom( aa ).replace('\n',''))
            
            fs.write_lines(fname, pdb_lines)
            self.parent().pymol.general_cmd(session, [('delete', 'all'), ('load', fname)])
            self.parent().pymol.remove_ligand_bonding(session)
            fs.safe_rm(fname)
            fname = str(Path(fname).with_suffix('.mol'))
            
            self.parent().pymol.general_cmd(session, [("save", fname)])
            
            if freeze_backbone:
                freeze_indices.extend( sr.backbone_indices(index=1))

            if freeze_hydrogens:
                for aa in sr.atoms:
                    if aa.element == 'H':
                        ha_parent = sr.closest_heavy_atom(aa)
                        dist = ha_parent.distance_to(aa)
                        constraint_lines.append(
                            f"   distance: {sr.get_atom_index(aa, 1)}, {sr.get_atom_index(ha_parent, 1)}, {dist:.3f}"         
                        )

            if constraints:
                for cst in constraints:
                    constraint_lines.extend( 
                        self.convert_constraint(sr, cst )
                    )

        else:

            if freeze_backbone:
                pass
    #            #TODO(CJ): add another constraint here
    #            pass
    #            input_lines.extend([
    #                '$fix',
    #               f'   atoms: {",".join(map(str,freeze_atoms))}',
    #                '$end'
    #            ])
    
            if constraints is not None:
                pass
    #            input_lines.append('$constrain')
    #            for cst in constraints:
    #                input_lines.append(
    #                f'   {cst[0]}: {", ".join(map(lambda ff: str(int(ff)), cst[-1]))}, {cst[1]:.2f}'
    #                )
    #            input_lines.append('$end')


            _parser.save_structure(fname, stru) 

        if freeze_indices:
            input_lines.extend([
                '$fix',
               f'   atoms: {",".join(map(str,sorted(freeze_indices)))}',
                '$end'
            ])


        if constraint_lines:
            input_lines.extend(["$constrain"] + constraint_lines + ["$end"])

        geom_input:str=None
        if input_lines:
            temp_path = Path(fname)
            geom_input = str(temp_path.parent / f"__geo_opt.inp")
            fs.write_lines(geom_input, input_lines)
            cmd_args.extend(["--input", geom_input ]) 

        cmd_args.extend([fname, "--opt" ])
        ext:str=fs.get_file_ext(fname)
        expected_output:str=f"xtbopt{ext}"
        results = self.env_manager_.run_command(
                self.config_.XTB_EXE,
                cmd_args)
        
        energy:float=None
        for ll in results.stdout.splitlines():
            if ll.find('TOTAL ENERGY') != -1:
                energy = float(ll.split()[3])

        session = self.parent().pymol.new_session()
        self.parent().pymol.general_cmd(session, [('set', 'retain_order', 1), ('load', expected_output)])
        df = self.parent().pymol.collect(session, 'memory', 'elem rank x y z'.split())
        self._remove_temp_files( str(Path(fname).parent) )
        self._remove_temp_files("./")
        for i, row in df.iterrows():
            point = (row.x, row.y, row.z)
            sr[i].coord = point
        
        return energy


    def _remove_temp_files(self, work_dir: str, extra_files:List[str]=None) -> None:
        """Removes the expected temp files created in the working directory of the input file. Deletes the file safely and 
        silently.

        Args:
            work_dir: Name of the directory to search for the files in.
            extra_files: A List[str] containing extra files to delete in the function.

        Returns:
            Nothing.
        """

        for fname in "charges wbo xtbrestart xtbtopo.mol xtbopt.log xtbopt.pdb xtbopt.mol xtbopt.xyz .xtboptok".split():
            fs.safe_rm(f"{work_dir}/{fname}")

        if extra_files:
            for ef in extra_files:
                fs.safe_rm(ef)

    def _check_valid_extension(self, fname: str) -> None:
        """Does the supplied file have a supported file extension? Logs an error
        and exits if not."""
        ext: str = fs.get_file_ext(str(fname))
        if ext in self.config_.SUPPORTED_EXTENSIONS:
            return

        _LOGGER.error(
            f"The supplied file '{fname}' has an unsupported extension of '{ext}'. Supported include {', '.join(self.config_.SUPPORTED_EXTENSIONS)}. Exiting..."
        )
        exit(1)


    @dispatch
    def convert_constraint( self, sr, cst: ResiduePairConstraint) -> List[str]:
    
        result = list()
        for (cst_name, child_cst) in cst.child_constraints:
            result.extend(self.convert_constraint(sr, child_cst ))

        return result

    @dispatch
    def convert_constraint( self, sr, cst:AngleConstraint) -> List[str]:
        #TODO(CJ): put in notes when the constraints are not present. Probably log a warning 
        result:List[str] = list()
        if not sr.has_atoms(cst.atoms):
            return list()
            print('error: TODO(CJ): make this better')
            exit( 0 )
        
        mapped_indices:List[int] = list(map(lambda aa: sr.get_atom_index(aa, indexing=1), cst.atoms))


        return [f"   angle: {mapped_indices[0]}, {mapped_indices[1]}, {mapped_indices[2]}, {cst.target_value:.3f}"]


    @dispatch
    def convert_constraint( self, sr, cst:DistanceConstraint) -> List[str]:
        result:List[str] = list()
        if not sr.has_atoms(cst.atoms):
            return list()
            print('error: TODO(CJ): make this better')
            exit( 0 )
        mapped_indices:List[int] = list(map(lambda aa: sr.get_atom_index(aa, indexing=1), cst.atoms))
        
        return [f"   distance: {mapped_indices[0]}, {mapped_indices[1]}, {cst.target_value:.3f}"]

