from typing import List
from abc import ABC

from enzy_htp import interface
from enzy_htp.geometry.sampling import equi_md_sampling
from enzy_htp import PDBParser, Structure

sp = PDBParser()
stru = sp.get_structure("xxx.pdb")
# general case
stru_esm: StructureEnsemble = equi_md_sampling(
    stru,
    engine="amber",
    equi_job_type="cpu",
    prod_time=100,
    cluster_job=True) # use kwarg to set each sampling

# constrain
constrain = [
    {"atoms" : [a1,a2],
     "constrain_value" : 2.0,}
]
stru_esm: StructureEnsemble = equi_md_sampling(
    stru,
    engine="amber",
    constrain=constrain,
    cluster_job=True)

# how to customize MD scheme?
parms = interface.amber.parameterize(stru)
step_1: MolDynStep = interface.amber.MolDynStep(stru, ntb=1, ntp=0,... , cluster_job=True)
step_2: MolDynStep = copy.deepcopy(step_1).set(ntb=2, ntp=1,... , cluster_job=True)
step_3: MolDynStep = copy.deepcopy(interface.amber.prod_step).set(stru, temp=2, cluster_job=True)
all_kinds_of_md_output = md_simulation(parms=parms, steps=[step_1, step_2, step_3], repeat=10) # could there be a better design; check engine consistency inside

class StructureEnsemble:
    def __init__(self, topology, top_format, coordinate, coord_format) -> None:
        # topology is a structure object with empty coord?
        # coordinate: List[List[Tuple[float, float, float]]]
        # both default stays in file?
        pass

    def get_stru_list(self) -> List[Structure]:
        pass

class MolDynStep(ABC): # abstract class defined in geometry? maybe in interface a file called requirements.py
    pass