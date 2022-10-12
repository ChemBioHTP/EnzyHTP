# type: ignore (for pylance)
'''
Scratch for what we want the final EnzyHTP API to be. (Use for design)

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2022-09-21
'''

# pylint: disable=all, E
import enzy_htp as eh

stru_parser = eh.structure.SeqParser()

def geno_enzy_htp_workflow(stru_parser: StructureParser):
    """
    science API planning for the GenoEnzyHTP application
    """
    candidate_seq = ['read_file_from_somewhere.fasta',]
    for seq in candidate_seq:
        raw_stru: Structure = stru_parser.get_structure(seq) # structural prediction seq -> stru
        
        stru_prepare_problems: List[StruProblem] = eh.detect_stru_problems(raw_stru, stru_parser) # pass stru_parser in this case problems are settled
        # 1
        prepped_stru = eh.prepare(raw_stru, stru_prepare_problems) # changes in the same stru but change name
        # 2
        preparer = eh.preparation.StructurePreparer(stru_prepare_problems) 
        preparer.add_problem(MissingLigand('path-to-ligand', information_of_docking, solver=eh.external_interface.Rosetta())) 
        preparer.delete_problem('stoichiometry') 
        preparer.problems['protonation'].set_ph(7.0)
        # will do the function version with **kwarg first





