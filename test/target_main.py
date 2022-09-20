# type: ignore (for pylance)
'''
Scratch for what we want the final EnzyHTP API to be. (Use for design)

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Date: 2022-08-30
'''

# pylint: disable=all, E
import enzy_htp as eh

stru_parser = eh.structure.structure_io.PDBParser()

def mutant_design_workflow(stru_parser: StructureParser):
    raw_stru = stru_parser.get_structure('a-path-to-file-or-a-object')

    #prepared_stru = prepare_structure(raw_stru)
    stru_prepare_problems: List[StruProblem] = eh.detect_stru_problems(raw_stru) # register preparer with required preparation steps (like a builder pattern)
    preparer = eh.preparation.StructurePreparer(stru_prepare_problems) # configure preparer to avoid massive function parameters ? should we use an object here? 
    preparer.add_problem(MissingLigand('path-to-ligand', information_of_docking, solver=eh.external_interface.Rosetta())) # the Problem classes have a common interface called solve that preparer do not need to know what kind of problem it is dealing. This design is because structural problem is easily expanding in type. In this design the problem carry the information of how to solve it in itself.
    preparer.delete_problem('stoichiometry') # one can edit the method of dealing with different problems
    preparer.problems['protonation'].set_ph(7.0) # currently preparer is nothing more than a special dictionary to contain problems

    prepared_stru = preparer.run(raw_stru) # this structure should be have complete topology data and initial coordinate

    muta_flag = eh.mutation.assign_mutation('pattern', prepared_stru) # this muta_flag is coupled with the specific struture as the mutator input
    mutator = eh.mutation.LeapMutator() # in the future there will be more mutators using different methods to mutate
    mutator.set_mutation(muta_flag)

    mutant_stru = mutator.run(prepared_stru) # its better to design unified interface for sciencific problem solver

    sampler = eh.sampling.MDEngine(eh.external_interface.Amber())
    sampler['production'].set_length(100, 'ns')
    
    stru_ensemble = sampler.run(mutant_stru)

    for stru in stru_ensemble:
        ele_stru_engine = eh.electronic.QMEngine(eh.external_interface.Gaussian())
        ele_stru_engine.software.scf_fail_solver.add_keyword('scf=qc') # these inteface objects can be configured with their specific problems
        ele_stru_engine.set_qm_region('byres. (resn. LIG around 5)')
        ele_stru_engine.set_qm_method('pbe1pbe', 'def2svp')
        ele_stru_data = e_engine.run(stru) # need to think more

def mutant_design_lambda_dynamic():
    '''
    planing about workflows that contain external processes that cover more than one htp level
    In this lambda dynamic case mutations are generated along the MD in the alchemical way.
    '''
    raw_stru = stru_parser.get_structure('a-path-to-file-or-a-object')
    stru_prepare_problems: List[StruProblem] = eh.detect_stru_problems(raw_stru)
    preparer = eh.preparation.StructurePreparer(stru_prepare_problems)
    prepared_stru = preparer.run(raw_stru)

    muta_flag = eh.mutation.assign_mutation('pattern', prepared_stru)
    alchemical_mutator = eh.mutation.alchemical.LambdaDynamicMutator() # alchemical mutator gives free energy differences of defined start/end state
    alchemical_mutator.set_mutation(muta_flag) # alchemical mutator also have mutator interfaces
    end_wt_stru = eh.structure.operation.undock(prepared_stru)
    
    mutant_stru, ddG_mut_start_end = alchemical_mutator.run(prepared_stru, end_wt_stru)

def mutant_design_workflow_w_object(parser, preparer, muta_flags, mutator, sampler, energy_engine, ranker):
    '''
    the use of objectized protocols decouples the workflow with its building blockings
    '''
    raw_stru = parser.get_struture('a-path-to-file-or-a-object')
    preparer.detect_stru_problems(raw_stru)
    prepared_stru = preparer.run(raw_stru)
    for muta_flag in muta_flags:
        mutant_stru = mutator.run(prepared_stru, muta_flag)
        stru_ensemble = sampler.run(mutant_stru)
        energy_features = energy_engine.run(stru_ensemble)
        ranker.add(muta_flag, energy_features)
    optimized_mutant = ranker.rank().get_best()

    return optimized_mutant

def mutant_ts_barrier_workflow():
    '''
    using free energy methods to calculate TS handles sampling and energy calculation the same time?
    '''
    pass

def mutant_design_workflow_wo_object():
    raw_stru = eh.structure.structure_io.parser_struture('a-path-to-file-or-a-object', 'pdb') # this way supporting a new format or changing an existing format will cause multiple code changes across layers
    problems = eh.preparation.detect_stru_problems(raw_stru)
    problems.append('docking') # code for configuring each step is coupled with the workflow code
    eh.preparation.prepare_stru(raw_stru, problems, docking_ligand_path='path-to-ligand', protonation_ph=7.0) # too many arguments are required. And it is coupled with the detect_stru_problem since it has to know which problem to provide information with
    pass

# a good idea in the avoid over engineering and use function as much as possible instead of classes. And those loosely defined interfaces can be made up by a well curated document.









