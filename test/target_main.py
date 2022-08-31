# type: ignore (for pylance)
'''
Scratch for what we want the final EnzyHTP API to be. (Use for design)

Author: shaoqz, <shaoqz@icloud.com>
Date: 2022-08-30
'''

# pylint: disable=all, E
import enzy_htp as eh

stru_parser = eh.structure.structure_io.PDBParser()

def mutant_design_workflow(stru_parser: StructureParser):
    raw_stru = stru_parser.get_structure('a-path-to-file-or-a-object')

    #prepared_stru = prepare_structure(raw_stru)
    stru_prepare_problems: List[StruProblem] = eh.detect_stru_problems(raw_stru) # register preparer with required preparation steps (like a builder pattern)
    preparer = eh.preparation.StructurePreparer(stru_prepare_problems) # configure preparer to avoid massive function parameters
    preparer.add_problem(MissingLigand('path-to-ligand', information_of_docking)) # the Problem classes have a common interface called solve that preparer do not need to know what kind of problem it is dealing. This design is because structural problem is easily expanding in type. In this design the problem carry the information of how to solve it in itself.
    preparer.delete_problem('stoichiometry') # one can edit the method of dealing with different problems
    preparer.problems['protonation'].set_ph(7.0) # currently preparer is nothing more than a special dictionary to contain problems

    prepared_stru = preparer.prepare(raw_stru) # this structure should be have complete topology data and initial coordinate












