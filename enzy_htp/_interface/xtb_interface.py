"""Defines an XTBInterface class that serves as a bridge for enzy_htp to utilize the xtb
semi-empirical quantum mechanics (QM) package. Uses the XTBConfig class found in enzy_htp/_config/xtb_config.py.
Supported operations include:
    + single point energy calculations
    + geometry optimization 


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-08-17
"""
from __future__ import annotations
import os
from pathlib import Path
from typing import Tuple, List, Dict
import numpy as np

from dataclasses import dataclass
from .base_interface import BaseInterface

from enzy_htp import _LOGGER
from enzy_htp.core import file_system as fs
from enzy_htp._config.xtb_config import XTBConfig, default_xtb_config
from enzy_htp.structure import Structure, Atom, StructureConstraint, PDBParser
from enzy_htp.structure.structure_region import (
    StructureRegion,
    create_region_from_residue_keys,
    create_region_from_full_stru
)
from enzy_htp.structure.structure_constraint import (
    CartesianFreeze,
    DistanceConstraint,
    AngleConstraint,
    DihedralConstraint,
    ResiduePairConstraint
)

from enzy_htp.electronic_structure import EletronicStructure

from enzy_htp.chemical import QMLevelOfTheory, MMLevelOfTheory, LevelOfTheory

from .handle_types import (
    QMSinglePointEngine,
    QMOptimizationEngine,
    QMResultEgg
)


from enzy_htp.core.job_manager import ClusterJob

from plum import dispatch
from enzy_htp import config as eh_config

@dataclass
class XTBQMResultEgg(QMResultEgg):
    """Class defining the ResultEgg for an XTB run."""
    charge_path:str
    wbo_path: str
    stru: Structure
    parent_job: ClusterJob


class XTBSinglePointEngine(QMSinglePointEngine):
    #TODO(CJ): documentation
    def __init__(self,
                interface,
                method: QMLevelOfTheory,
                region: StructureRegion,
                keep_geom: bool,
                name: str,
                cluster_job_config: Dict,
                keep_in_file: bool,
                work_dir: str,
                ):
        self._parent_interface = interface
        self._method = method
        self._region = region
        self._name = name
        self._cluster_job_config = cluster_job_config
        self._keep_in_file = keep_in_file
        self._work_dir = work_dir
        self._geo_opt = False
        self._constraints = list()

    @property
    def method(self) -> QMLevelOfTheory:
        """Getter for the QMLevelOfTheory dictating the Engine's theoretical behavior."""
        return self._method

    @property
    def geo_opt(self) -> bool:
        """Is this a QMGeometryOptimizationEngine?"""
        return self._geo_opt

    @property
    def parent_interface(self) -> XTBInterface:
        """Getter for the Engine's parent interface."""
        return self._parent_interface

    @property
    def work_dir(self) -> str:
        return self._work_dir

    @property
    def region(self) -> StructureRegion:
        """Getter for the StructureRegion """
        return self._region

    @property
    def engine(self) -> str:
        """Name of the engine. Hardcoded to 'xtb'"""
        return "xtb"

    @property
    def constraints(self) -> List[StructureConstraint]:
        """Getter for the StructureConstraints owned by the Engine."""
        return self._constraints

    def make_job(self, stru: Structure) -> Tuple[ClusterJob, XTBQMResultEgg]:
        """The method that makes a ClusterJob that runs the QM."""
        # 1. input
        if not isinstance(stru, Structure):
            _LOGGER.error("The supplied variable stru MUST be a Structure()")
            raise TypeError
        assert False, "TODO(CJ): translate for xtb"
        # 2. make .gjf file
        fs.safe_mkdir(self.work_dir)
        #temp_gjf_file, gchk_path = self._make_gjf_file(stru)

        # 3. make cmd
        spe_cmd, gout_path = self.parent_interface.make_gaussian_cmd(temp_gjf_file)

        # 4. assemble ClusterJob
        cluster = self.cluster_job_config["cluster"]
        res_keywords = self.cluster_job_config["res_keywords"]
        env_settings = cluster.G16_ENV["CPU"]
        sub_script_path = fs.get_valid_temp_name(f"{self.work_dir}/submit_{self.name}.cmd")
        job = ClusterJob.config_job(
            commands = spe_cmd,
            cluster = cluster,
            env_settings = env_settings,
            res_keywords = res_keywords,
            sub_dir = "./", # because path are relative
            sub_script_path = sub_script_path
        )
        job.mimo = { # only used for translate clean up
            "temp_gin": [temp_gjf_file],
        }

        # 5. make result egg
#        result_egg = GaussianQMResultEgg(
#            gout_path = gout_path,
#            gchk_path = gchk_path,
#            stru=stru,
#            parent_job = job,
#        )

        return (job, result_egg)

        assert False

    def run(self, stru: Structure) -> EletronicStructure:
        """Method that actually runs the XTB single point calculation. Returns results of calculation as an ElectronicStructure() object."""
        if not isinstance(stru, Structure):
            _LOGGER.error("Supplied stru is not of type Structure()")
            raise TypeError()
        
        fs.safe_mkdir(self.work_dir)

        sr = self.region
        if sr is None:
            assert False, "Need to fix this logic"
        
        run_info:Dict = self.parent_interface.setup_xtb_run(
            sr,
            constraints=self.constraints,
            lot=self.method,
            geo_opt=self.geo_opt,
            work_dir=self.work_dir,
        )

        start_dir:str = os.getcwd()
        os.chdir(run_info['work_dir'])

        results = self.parent_interface.env_manager_.run_command(
                self.parent_interface.config()['XTB_EXE'],
                run_info['args'])

        os.chdir(start_dir)

        lines:List[str] = fs.lines_from_file(run_info['xtb_outfile'])

        return  EletronicStructure(
            energy_0 = self.parent_interface.parse_spe(lines),
            geometry = stru,
            mo = '',
            mo_parser = '',
            source="xtb",
        )

    def translate(self, result_egg: XTBQMResultEgg) -> EletronicStructure:
        """TODO(CJ)"""

        assert False

class XTBOptimizationEngine(XTBSinglePointEngine, QMOptimizationEngine):
    def __init__(self,
                interface,
                method: QMLevelOfTheory,
                region: StructureRegion,
                constraints: List[StructureConstraint],
                keep_geom: bool,
                name: str,
                cluster_job_config: Dict,
                keep_in_file: bool,
                work_dir: str,
                ):
        
        self._constraints = constraints
        self._parent_interface = interface
        self._method = method
        self._region = region
        self._name = name
        self._cluster_job_config = cluster_job_config
        self._keep_in_file = keep_in_file
        self._work_dir = work_dir
        self._geo_opt = True 


    def run(self, stru: Structure) -> EletronicStructure:
        """Method that actually runs the XTB single point calculation. Returns results of calculation as an ElectronicStructure() object."""
        if not isinstance(stru, Structure):
            _LOGGER.error("Supplied stru is not of type Structure()")
            raise TypeError()
        
        fs.safe_mkdir(self.work_dir)

        sr = self.region
        if sr is None:
            sr = create_region_from_full_stru( stru )
        
        run_info:Dict = self.parent_interface.setup_xtb_run(
            sr,
            constraints=self.constraints,
            lot=self.method,
            geo_opt=self.geo_opt,
            work_dir=self.work_dir,
        )

        start_dir:str = os.getcwd()
        os.chdir(run_info['work_dir'])

        results = self.parent_interface.env_manager_.run_command(
                self.parent_interface.config()['XTB_EXE'],
                run_info['args'])

        os.chdir(start_dir)

        self.parent_interface.update_coords(sr, run_info['expected_outfile'])

        lines:List[str] = fs.lines_from_file(run_info['xtb_outfile'])

        return  EletronicStructure(
            energy_0 = self.parent_interface.parse_spe(lines),
            geometry = stru,
            mo = '',
            mo_parser = '',
            source="xtb",
        )



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

    def setup_xtb_run(self,
        stru:Union[StructureRegion],
        lot:QMLevelOfTheory, 
        geo_opt:bool = False,
        charge:int = None,
        spin:int = 1,
        constraints:List[StructureConstraint]=None,
        n_iter: int = None,
        n_proc: int = None,
        work_dir:str=None
        ) -> Dict:
        """This method does the majority of the leg work before an xtb calculation is performed. All relevant variables are stored in 
        the result dict(). Note that the same method is used for both single point energy calculations and geometry optimizations. 

        Args:
            stru:
            lot:
            geo_opt:
            charge:
            spin:
            constraints:
            n_iter:
            n_proc:
            work_dir:

        Returns:

        """

        if not isinstance(stru, StructureRegion):
            _LOGGER.error("The supplied geometry must be of type StructureRegion()")
            raise TypeError()
    
        self.check_valid_lot(lot)

        if charge is None:
            charge = stru.get_net_charge()

        if n_iter is None:
            n_iter = self.config().N_ITER

        if n_proc is None:
            n_proc = self.config().N_PROC

        result = dict()

        coord_file:str = self.write_coord_file(stru, work_dir)
        coord_path = Path(coord_file)
        result['coord_file'] = coord_file
        
        work_dir_abs:str=f"{Path(work_dir).absolute()}"
        result['work_dir'] = work_dir_abs 
        
        xtb_outfile:str=f"{work_dir_abs}/xtb_log.txt"
        result['xtb_outfile'] = xtb_outfile

        args:List[str] = ["--chrg", str(charge), "--iterations",
                        str(n_iter), "--parallel",
                        str(n_proc), "--norestart", ]
        if geo_opt:
            if constraints:
                xtb_inp_file:str=f"{work_dir_abs}/xtb_settings.inp"
                self.write_inp_file(xtb_inp_file, stru, constraints)
                result['xtb_inp_file'] = xtb_inp_file
            
                args.extend(["--input", xtb_inp_file])

            args.extend(["--opt", "normal"]) #TODO(CJ): fix this later

            result['expected_outfile'] = str( coord_path.parent / f"xtbopt{coord_path.suffix}" )

        else:
            args.append( "--sp" )

        args.extend([str(coord_path.absolute()), ">", xtb_outfile])
       
        result['args'] = args

        return result

    def parse_spe(self, lines:List[str]) -> float:
        """Given the output lines from an xtb run, parse and return the single 
        point energy as a float(). Throws an error if not found.

        Args:
            lines: The List[str] with the results of the xtb run.

        Returns:
            The single point energy as a float in hartree units. 

        """
        for ll in lines:
            if ll.find('TOTAL ENERGY') != -1:
                return float(ll.split()[3])
        else:
            raise TypeError("No total energy value in xtb output!")


    def write_coord_file(self, sr:StructureRegion, work_dir:str) -> str:
        """Given a StructureRegion and work_dir, write the xtb coords of the region. The path it will be written
        to is work_dir/xtb_coords.mol.

        Args:
            sr: The StructureRegion describing the QM region of the system.
            work_dir: The str() of the directory to write the coords file.

        Returns:
            The path of the coords .mol file.

        """
        temp_pdb:str = f"{work_dir}/temp_xtb_pdb.pdb"
        coord_file:str = f"{work_dir}/xtb_coords.mol"
        _parser = PDBParser()

        pdb_lines:List[str] = list()
        for aidx,aa in enumerate(sr.atoms):
            aa.idx = aidx
            pdb_lines.append( _parser._write_pdb_atom( aa ).replace('\n',''))
        
        fs.write_lines(temp_pdb, pdb_lines)
        session = self.parent().pymol.new_session()
        self.parent().pymol.general_cmd(session, [('delete', 'all'), ('load', temp_pdb)])
        self.parent().pymol.remove_ligand_bonding(session)
        fs.safe_rm(temp_pdb)
        
        self.parent().pymol.general_cmd(session, [("save", coord_file)])

        return coord_file


    def build_single_point_engine(
        self,
        # calculation config
        method: QMLevelOfTheory = "default",
        region: StructureRegion = None,
        keep_geom:bool = "default",
        # calculation config (alternative)
        name: str = "default",
        cluster_job_config: Dict = "default",
        keep_in_file: bool = False,
        work_dir: str = None,) -> XTBSinglePointEngine:
        """Simple constructor that builds the XTBSinglePointEngine. Designed to work with quantum.single_point() API."""

        return XTBSinglePointEngine(
            interface=self,
            method=method,
            region=region,
            keep_geom=keep_geom,
            name=name,
            cluster_job_config=cluster_job_config,
            keep_in_file=keep_in_file,
            work_dir=work_dir
        )


    def build_optimization_engine(
        self,
        method:QMLevelOfTheory,
        region:StructureRegion,
        constraints:List[StructureConstraint] = None,
        keep_geom: bool = "default", 
        name: str = "default",
        cluster_job_config: Dict = "default",
        keep_in_file: bool = False,
        work_dir:str = None,) -> XTBOptimizationEngine:
        """Simple constructor that builds the XTBOptimizationEngine. Designed to work with the quantum.optimization() API."""
        
        return XTBOptimizationEngine(
            interface=self,
            method=method,
            region=region,
            constraints=constraints,
            keep_geom=keep_geom,
            name=name,
            cluster_job_config=cluster_job_config,
            keep_in_file=keep_in_file,
            work_dir=work_dir
        )
            


    def write_inp_file(self, out_path:str,
                            stru: StructureRegion,
                            constraints:List[StructureConstraint],
                            ) -> str:
        """Writes a .inp file

        Args:

        Returns:

        """
        constraint_lines:List[str] = list()
        frozen_indices:List[int] = list()

        for cst in constraints:
            if cst.is_cartesian_freeze():
                for aa in cst.atoms:
                    if stru.has_atom(aa):
                        frozen_indices.append( stru.get_atom_index(aa, indexing=1) )

            else:
                constraint_lines.extend( 
                    self.convert_constraint(stru, cst )
                )
        
        inp_lines = list()
        if frozen_indices:
            frozen_indices.sort()
            fi_str = ",".join(map(str,frozen_indices))
            inp_lines.extend([            
            "$fix",
            f"   atoms: {fi_str}",
            "$end"
            ])

        if constraint_lines:
            inp_lines.extend(["$constrain"] + constraint_lines + ["$end"])

        fs.write_lines(out_path, inp_lines )

        return out_path

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

    def check_valid_lot(self, lot: QMLevelOfTheory) -> None:
        """Helper method that checks if the given QMLevelOfTheory is valid. Throws an error if not.
        Args:
            lot: The QMLevelOfTheory to check.

        Returns:
            Nothing.
        """
        if lot.method and lot.method not in self.config()['SUPPORTED_XTB_THEORY_LEVELS']:
            _LOGGER.error(f"The method {lot.method} is not supported. Allowed values are: {', '.join(self.config()['SUPPORTED_XTB_THEORY_LEVELS'])}")
            raise TypeError()

        if lot.solv_method and lot.solv_method not in self.config()['SOLVATION_METHODS']:
            _LOGGER.error(f"The solvation method {lot.solv_method} is not supported. Allowed values are: {', '.join(self.config()['SOLVATION_METHODS'])}")
            raise TypeError()

        if lot.solv_method == 'ALPB':
            if lot.solvent not in self.config()['ALPB_SOLVENTS']:
                _LOGGER.error(f"The solvent {lot.solvent} is not supported in the ALPB model. Allowed values are: {', '.join(self.config()['ALPB_SOLVENTS'])}")
                raise TypeError()
        
        elif lot.solv_method == 'GBSA':
            if lot.solvent not in self.config()['ALPB_SOLVENTS']:
                _LOGGER.error(f"The solvent {lot.solvent} is not supported in the GBSA model. Allowed values are: {', '.join(self.config()['ALPB_SOLVENTS'])}")
                raise TypeError()


    def update_coords(self, sr, coord_file:str) -> None:
        """TODO(CJ)"""
        session = self.parent().pymol.new_session()
        df:pd.DataFrame = self.parent().pymol.collect(session, coord_file, "x y z rank elem".split())
        df.sort_values(by='rank', inplace=True)
        #TODO(CJ): add some checks in here
        for aa, (i, row) in zip(sr.atoms, df.iterrows()):
            assert aa.element == row['elem']
            aa.coord = np.array([row['x'], row['y'], row['z']])

xtb_interface = XTBInterface(None, eh_config._xtb)
"""The singleton of XTBInterface() that handles all XTB related operations in EnzyHTP
Instantiated here so that other _interface subpackages can use it.
An example of this concept this AmberInterface used Gaussian for calculating the RESP charge
so it imports gaussian_interface from here."""
