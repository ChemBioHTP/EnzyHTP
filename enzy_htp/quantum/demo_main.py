"""a design sketch of this module"""
from typing import List, Tuple
from abc import ABC

from enzy_htp import interface
from enzy_htp.geometry.sampling import equi_md_sampling
from enzy_htp.quantum import single_point
from enzy_htp.structure import StructureEnsemble
from enzy_htp import PDBParser, Structure


class EletronicStructure:
    # attribute
    energy_0: float
    geometry: Structure
    mo: str
    mo_parser: callable
    # deduced
    mo_coeff: List[List[float]]
    mo_occ: List[float]
    basis_set: str

class LevelofTheory:
    _type: str

class QMLevelofTheory(LevelofTheory):
    basis_set: str
    method: str
    solvent: str
    solv_method: str

class MMLevelofTheory(LevelofTheory):
    force_field: str
    ligand_method: str

sp = PDBParser()
stru = sp.get_structure("xxx.pdb")
param_method = interface.amber.build_md_parameterizer()
stru_esm: List[StructureEnsemble] = equi_md_sampling(
                                        stru,
                                        param_method,
                                        cluster_job_config={
                                            "account" : "csb_gpu_acc",
                                            "partition" : "turing"})

# single point can take 1 stru
ele_stru: EletronicStructure = single_point(
    stru,
    QMLevelofTheory(
        basis_set="def2-svp", method="b3lyp-d3",
        solvent="h2o", solv_method="smd",
    ),
    use_symmetry=False,)

# single point can take an ensemble
ele_stru: EletronicStructure = single_point(
    stru_esm,
    QMLevelofTheory(
        basis_set="def2-svp", method="b3lyp-d3",
        solvent="h2o", solv_method="smd",
    ),
    use_symmetry=False,)

# single point can take an ensemble and calculate only a region of it.
ele_stru: EletronicStructure = single_point(
    stru_esm,
    regions=["resi 101+254"],
    region_methods = [
        QMLevelofTheory(
            basis_set="def2-svp", method="b3lyp-d3",
            solvent="h2o", solv_method="smd",
        ),
    ],
    engine="g16",)

# single point can take an ensemble and calculate in a multiscale manner.
ele_stru: EletronicStructure = single_point(
    stru_esm,
    regions=["resi 101+254", "else"],
    region_methods = [
        QMLevelofTheory(
            basis_set="def2-svp", method="b3lyp-d3",
            solvent="h2o", solv_method="smd",
        ),
        MMLevelofTheory(
            force_field=["ff14sb", "gaff2"],
            ligand_method="gaff2-resp"
        ),
    ],
    embedding_method="mechanical",
    engine="chemshell",) # electrostatic, solid-state, ...


