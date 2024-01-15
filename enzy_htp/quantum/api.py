"""Define functions for quantum mechanics calculation of Structure(). These functions derives electronic
structure based on atomic coordinates using DFT or other QM methods.
The electronic structure can then be used to calculate energy, forces, and other useful descriptors.

Science API:
    + single_point()
    + optimize()

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-10-28"""
from typing import Union, List, Dict, Callable

from enzy_htp import interface
from enzy_htp._interface.handle_types import QMSinglePointEngine, QMOptimizationEngine
from enzy_htp.electronic_structure import ElectronicStructure
from enzy_htp.structure import (
    Structure,
    StructureEnsemble,
    create_region_from_selection_pattern
)
from enzy_htp.structure.structure_operation import init_charge
from enzy_htp.chemical import QMLevelOfTheory, LevelOfTheory
from enzy_htp.core.logger import _LOGGER
import enzy_htp.core.job_manager as armer

def single_point(
        stru: Union[Structure, StructureEnsemble],
        engine: str, # always acknowledge the engine
        # single region case option
        method: QMLevelOfTheory = None, # single region has to be QM
        # multi region case option
        regions: List[str]= None,
        region_methods: List[LevelOfTheory]= None,
        capping_method: str = "res_ter_cap",
        embedding_method: str= "mechanical", # TODO probably not a good default choice
        parallel_method: str="cluster_job",
        cluster_job_config: Dict= None,
        job_check_period: int= 210, # s
        job_array_size: int= 20,
        work_dir: str="./QM_SPE",
        keep_in_file: bool=False,
        ) -> List[ElectronicStructure]:
    """The QM single point calculation. This function calculates the molecular orbitals (MOs)
    give a molecule with a specific geometry. If an ensemble of geometry is give, 
    calculation is applied to each snapshot in this ensemble.
    Args:
        stru:
            the target molecule of the calculation represented as Structure()
            It can also be an ensemble of structures as StructureEnsemble()
            and in this case, each geometry in this ensemble will be calculated.
        engine:
            the QM or QM/MM engine.
        method:
            the level of theory of this calculation as a LevelOfTheory().
            This is used when there is only 1 or less region specified.
        regions:
            This option allows you to define different region and apply different
            level of theory to each region.
            e.g.: defining 2 regions and perform QM/MM.
        region_methods:
            The level of theory of each region.
            This is used when more than 1 region is specified.
            The region and the method is align based on the order.
        capping_method:
            the free valence capping method.
        embedding_method:
            The embedding method of multiscale simulation.
            This is used when more than 1 region is specified.
        parallel_method:
            the method to parallelize the multiple runs when more
            than 1 geometry is in the input StructureEnsemble
        job_check_period:
            the time cycle for update job state change (Unit: s)
        job_array_size:
            how many jobs are allowed to submit simultaneously. (0 means all -> len(inp))
            (e.g. 5 for 100 jobs means run 20 groups. All groups will be submitted and
            in each group, submit the next job only after the previous one finishes.)
        work_dir:
            the working dir that contains all the files in the SPE process
        keep_in_file:
            whether keep the input file of the calculation
    Returns:
        A list of ElectronicStructure for each frame
        [ElectronicStructure(), ...]
    
    Details:
        The quantum mechanics models molecules by solving their corresponding
    Schrodinger equation. This will yield a relationship between the electronic
    motion/state (representated as wavefunctions or orbitals) and energy. c This enables geometry optimization and in a geometry with well-defined
    meaning, predicting the molecule's response to energy related events. Besides the
    energy itself, the ground state wavefunction is also useful to derivate all desired
    properities and assess how the given geometry will respond to a certain field. For
    example, the dipole moment of the breaking bond can be calculated to be a useful
    descriptor for determine how much the internel electric field of the enzyme stablize
    the reaction.
        The level of theory defines the method used to solve the Schrodinger equation and
    basis function set that used for the description of the MOs. An inappropriate selection
    of level of theory may leads to completely wrong result. The selection of the level
    of theory for different tasks is widely discussed/benchmarked in the internet. Search
    for benchmark papers or selections of systems close to your need. Here are several good
    summary (in Chinese): http://sobereva.com/336, http://sobereva.com/272.
        It is also possible to apply different level of theory in different regions of a
    structure. For example, use DFT for the active site, and use MM for the rest of the
    enzyme. In this case, another problem emergies, that is, how to consider the influence
    from one region to the other, the embedding probem.

    Engines that we plan (or already) to support:
        Gaussian:
            + best optimization method, fast convergence.
            + stable and robust
            + huge community
        XTB:
            + gfn-xTB series methods fast and relatively accurate.
            + crest conformational search tool.
        Orca:
            + DLPNO-CCSD(T) practicable gold standard method
            + NEB for double ended TS search
        pySCF:
            + python interface
            + fastest SCF
        QChem:
            + EDA
        ChemShell:
            + Multiscale QM/MM
            + python interface"""
    
    # init
    # dispatch: stru -> stru_esm
    if isinstance(stru, Structure):
        stru_esm = StructureEnsemble.from_single_stru(stru)
    elif isinstance(stru, StructureEnsemble):
        stru_esm = stru
    else:
        _LOGGER.error(f"only accept Structure() or StructureEnsemble(). Got: {stru}")
        raise TypeError

    # dispatch: region, region_methods, method -> qm_engine
    if regions is None:
        regions = []
    if region_methods is None:
        region_methods = []

    if len(regions) != len(region_methods):
        if not len(regions) == 1:
            _LOGGER.error(f"Number of region methods ({len(region_methods)}) do not match number of regions ({len(regions)}).")
            raise ValueError
    
    if len(regions) < 2:
        # QM only
        # init
        # qm_engine
        qm_engine_ctor = SINGLE_REGION_SINGLE_POINT_ENGINE[engine]
        # qm_method
        if region_methods:
            if not isinstance(region_methods[0], QMLevelOfTheory):
                _LOGGER.error(f"Only 1 or less region specified. Have to be a QMLevelOfTheory. Got: {region_methods[0]}")
                raise TypeError
            qm_method = region_methods[0]
        elif not isinstance(method, QMLevelOfTheory):
            _LOGGER.error(f"Have to specify a QM level of theory in `method` or `region_methods` Got: {method}")
            raise TypeError
        else:
            qm_method = method
        # stru_region
        if not regions:
            # whole
            qm_region = None
        else:
            qm_region = create_region_from_selection_pattern(
                stru_esm.topology,
                regions[0],
                capping_method)
            init_charge(qm_region)
        
        qm_engine: QMSinglePointEngine = qm_engine_ctor(
                                            region=qm_region,
                                            method=qm_method,
                                            keep_geom=True,
                                            cluster_job_config=cluster_job_config,
                                            work_dir=work_dir,
                                            keep_in_file=keep_in_file,
                                        )
    else:
        # multiscale
        raise Exception("TODO")
        qm_engine = MULTI_REGION_SINGLE_POINT_ENGINE[engine]
        regions = [create_region_from_selection_pattern(stru_esm.topology, i) for i in regions]
        qm_engine = qm_engine_ctor(
                        regions=regions,
                        region_methods=region_methods,
                        embedding_method=embedding_method,
                        keep_geom=True,
                        cluster_job_config=cluster_job_config,
                        work_dir=work_dir,
                        keep_in_file=keep_in_file,
                    )

    # parallel methods
    if parallel_method == "cluster_job":
        if not cluster_job_config:
            _LOGGER.error("cluster_job is used but cluster_job_config is not given! "
                        "You need to at least specify the account and partition. "
                        "See test/geometry/test_sampling.py::test_equi_md_sampling_lv1() for an example.")
            raise ValueError
        result = _parallelize_qm_with_cluster_job(
            stru_esm, qm_engine,
            job_check_period,
            job_array_size,
        )
    elif parallel_method is None:
        result = _serial_qm(
            stru_esm, qm_engine,
        )
    else:
        _LOGGER.error(f"{parallel_method} is not a supported parallel_method")
        raise ValueError        

    return result


def optimize(stru: Structure):
    raise Exception("TODO")


SINGLE_REGION_SINGLE_POINT_ENGINE: Dict[str, Callable] = {
    "g16": interface.gaussian.build_single_point_engine,
    "gaussian": interface.gaussian.build_single_point_engine,
    "gaussian16": interface.gaussian.build_single_point_engine
}
"""constructors of QMSinglePointEngine for single_point() when only
1 region is used. full QM."""


MULTI_REGION_SINGLE_POINT_ENGINE: Dict[str, Callable] = {
    "chemshell": "TODO",
    "chemsh": "TODO",
    "g16": interface.gaussian.build_qmmm_single_point_engine,
    "gaussian": interface.gaussian.build_qmmm_single_point_engine,
    "gaussian16": interface.gaussian.build_qmmm_single_point_engine,
}
"""constructors of QMSinglePointEngine for single_point() when more than
1 region is used. multiscale QM."""


def _parallelize_qm_with_cluster_job(
        stru_esm: StructureEnsemble,
        qm_engine: Union[QMSinglePointEngine, QMOptimizationEngine],
        job_check_period: int,
        array_size: int,
        ) -> List[ElectronicStructure]:
    """The QM parallelization method: cluster_job.
    This method will utilize ARMer@EnzyHTP and make each QM calculation a ClusterJob and
    parallalize them in a job array
    
    this potentially should work for both single point and opt.
    modify when add finish opt function. In that case, need to update type hinting"""
    result = []
    qm_job_list = []
    qm_result_eggs = []
    # 1. prep jobs
    for stru in stru_esm:
        qm_job, output = qm_engine.make_job(stru)
        qm_job_list.append(qm_job)
        qm_result_eggs.append((qm_engine, output))
    
    # 2. wait to end
    armer.ClusterJob.wait_to_array_end(
        qm_job_list,
        period=job_check_period,
        array_size=array_size,)
    
    # 3. translate eggs
    for engine, egg in qm_result_eggs:
        result.append(engine.translate(egg))

    return result


def _serial_qm(
        stru_esm: StructureEnsemble,
        qm_engine: Union[QMSinglePointEngine, QMOptimizationEngine],
        ) -> List[ElectronicStructure]:
    """The QM serial running method
    This method runs QMs in a serial manner locally."""
    result = []
    # 1. run jobs
    for stru in stru_esm:
        output = qm_engine.run(stru)
        result.append(output) 
    
    return result
