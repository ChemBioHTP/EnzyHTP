"""Driver for the reactive docking functionality available in enzy_htp. The only function that that should be 
called is dock_reactants(). All others are implementation functions that SHOULD NOT be used. 

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-07-28
"""
import os
from pathlib import Path
from collections import defaultdict
from typing import List, Tuple, Set, Dict
from copy import deepcopy

import numpy as np
import pandas as pd

from enzy_htp import interface, config, _LOGGER
import enzy_htp.structure.structure_operation as stru_oper
from enzy_htp.structure.structure_constraint import StructureConstraint, CartesianFreeze

import enzy_htp.chemical as chem
from enzy_htp.structure import PDBParser, Mol2Parser, Structure, Ligand, translate_structure, Atom
from enzy_htp.core import file_system as fs
from enzy_htp.quantum import single_point
from enzy_htp.quantum import optimize as qm_optimize 

def dock_reactants(structure: Structure,
                   ligands: List[Ligand],
                   constraints: List[StructureConstraint] = None,
                   n_struct: int = 100,
                   cst_energy: float = None,
                   use_qm: bool = True,
                   clash_cutoff: int = 3,
                   max_sasa_ratio:float = 0.60,
                   cluster_distance: float = 2.0,
                   contact_threshold:float=3.0,
                   box_size:float=50.0,
                   move_distance:float=40.0,
                   transform_angle:float=360,
                   transform_cycles:int=1000,
                   transform_repeats:int=3,
                   transform_temperature:int=5,
                   fr_repeats:int=1,
                   grid_width:float=50.0,
                   rng_seed: int = 1996,
                   work_dir: str = None,
                   save_work_dir: bool = True,
                   ) -> None:
    """Takes a Structure() containing Ligand() objects and tries to create a complex with an optimized geometry 
    consistent with that described in the supplied constraints. Does all work inplace on the supplied Structure().
    Below is the workflow of the function:

    1. Generate constrained geometries with RosettaLigand()
    2. Filter down and select a target geometry
        a. Filter out if too many clashes
        b. Filter out if ligands have too much SASA
        c. Filter out if ligands do not satisfy constraints
        d. Select geometry with:
            i. Lowest QM energy, if use_qm is True
            ii. Lowest RosettaLigand energy, if use_qm is False
    3. Rosetta FastRelax minimization 1
        a. use constraints
        b. allow backbone flexibility in active site
    4. constrained QM minimization of active site with xtb (if use_qm is True)
    
    Args:
        structure: The Structure() to perform reactive docking on. Contains both Ligand()'s and protein chains.
        constraints: A List[StructureConstraint] describing the reaction geometry.
        n_struct: How many geometries should we make with RosettaLigand? Default is 100.
        cst_energy: The energy penalty cutoff for keeping a given geometry. In REU and assuming Rosetta constraint scoring.
        use_qm: Shoud QM be used for both geometry filtration and geometry minimization? Default is True.
        freeze_alphafill: If a Ligand() was placed with AlphaFill, should it be frozen during initial geometry sampling. Default is True.
        clash_cutoff: How many heavy atom clashes are allowed in a given geometry? Default is 3.
        max_sasa_ratio: What percentage of a Ligand()'s surface area can be solvent accessible?  Default is 0.60.
        cluster_distance: How close should a residue be to the Ligand()'s to be included in the QM region? Default is 2.0 Angstroms.
        rng_seed: Random number generation integer seed. Default value is 1996.
        work_dir: Where is the work being done/where should the scratch files be made?
        save_work_dir: Should temporary files be saved? Default is True.

    Returns:
        Nothing.
    """

    if work_dir is None:
        work_dir = config["system.SCRATCH_DIR"]

    if cst_energy is None:
        cst_energy = len(constraints)*2000.0

    fs.safe_mkdir(work_dir)

    translate_structure(structure, end_naming='rosetta')
    
    param_files:List[str] = parameterize_structure(structure, work_dir)

    for ligand in ligands:
        if id(ligand.root()) != id(structure):
            err_msg:str=f"The supplied ligand {ligand} is not a child of the supplied structure!"
            raise TypeError(err_msg)
    

    ligands.reverse()
    while ligands:
        ligand = ligands.pop()
        relevant_csts:List[StructureConstraint] = list()
        for cst in constraints:
            if cst.is_constraining(ligand):
                for ll in ligands:
                    if cst.is_constraining(ll):
                        break
                else:
                    relevant_csts.append(cst)


        dock_ligand(structure,
                        ligand,
                        relevant_csts,
                        param_files,
                        n_struct,
                        clash_cutoff,
                        max_sasa_ratio,
                        cst_energy,
                        cluster_distance,
                        contact_threshold,
                        box_size,
                        move_distance,
                        transform_angle,
                        transform_cycles,
                        transform_repeats,
                        transform_temperature,
                        1,#fr_repeats,
                        grid_width,
                        use_qm,
                        rng_seed,
                        f"{work_dir}/rosetta_docking/")


    sp = PDBParser()

    mm_minimization(structure, constraints, param_files, contact_threshold, False, fr_repeats, rng_seed, f"{work_dir}/rosetta_minimization/")

    if use_qm:
        qm_minimization(structure, constraints, cluster_distance, False, work_dir)


    mm_minimization(structure, constraints, param_files, contact_threshold, True, fr_repeats, rng_seed, f"{work_dir}/rosetta_minimization/")


    if use_qm:
        qm_minimization(structure, [], cluster_distance, True, work_dir)

    translate_structure(structure, start_naming='rosetta')

    if not save_work_dir:
        _LOGGER.info(f"save_work_dir set to False! Deleting {work_dir}")
        fs.safe_rmdir( work_dir )


def mm_minimization(structure:Structure,
                constraints:List[StructureConstraint],
                param_files:List[str],
                contact_threshold:float,
                ramp_constraints:bool,
                fr_repeats:int,
                rng_seed:int,
                work_dir:str) -> None:
    """
    """

    fs.safe_mkdir(work_dir)
    xml_script:str = create_minimization_xml(work_dir)
    pdb_file:str=f"{work_dir}/start.pdb"
    _sp = PDBParser()
    _sp.save_structure(pdb_file, structure)
    cst_file:str = interface.rosetta.write_constraint_file(structure, constraints, work_dir) #TODO(CJ): look at this; wrong constraint types!!

    sele = list()
    for cst in constraints:
        for atom in cst.atoms:
            parent_residue = atom.parent
            sele.append(f"{parent_residue.idx}{parent_residue.parent.name}")

    sele = ",".join(list(set(sele)))
        

    #TODO(CJ): make the sele here
    options_file:str = create_minimization_options(
        work_dir,
        pdb_file,
        xml_script,
        param_files,
        rng_seed,
        sele,
        contact_threshold,
        ramp_constraints,
        cst_file,
        fr_repeats
    )

    start_dir: str = os.getcwd()
    
    os.chdir(work_dir)
    
    try:
        interface.rosetta.run_rosetta_scripts([f"@{Path(options_file).name}"])
    except:
        pass


    os.chdir(start_dir)

    scores_file:str = f"{work_dir}/minimization_scores.sc"
    
    fs.check_file_exists(scores_file, exit_script = False )

    df: pd.DataFrame = interface.rosetta.parse_score_file(scores_file)

    df['description'] = df.apply(lambda row: f"{work_dir}/{row.description}.pdb", axis=1)
    
    energy_key='total_score'

    infile:str=df.sort_values(by=energy_key).description.to_list()[0]

    _parser = PDBParser()
    ref_stru = _parser.get_structure(infile)
    stru_oper.update_residues(structure, ref_stru)

    for cst in constraints:
        cst.change_topology(structure)

def dock_ligand(structure:Structure,
                    ligand:Ligand,
                    constraints:List[StructureConstraint],
                    param_files:List[str],
                    n_struct:int,
                    clash_cutoff:float,
                    max_sasa_ratio:float,
                    cst_energy:float,
                    cluster_distance:float,
                    contact_threshold:float,
                    box_size:float,
                    move_distance:float,
                    transform_angle:float,
                    transform_cycles:int,
                    transform_repeats:int,
                    transform_temperature:float,
                    fr_repeats:int,
                    grid_width:float,
                    use_qm:bool,
                    rng_seed:int,
                    work_dir:str) -> None:


    fs.safe_mkdir(work_dir)
    xml_script:str = create_docking_xml(work_dir)
    pdb_file:str=f"{work_dir}/start.pdb"
    _sp = PDBParser()
    _sp.save_structure(pdb_file, structure)

    cst_file:str = interface.rosetta.write_constraint_file(structure, constraints, work_dir) #TODO(CJ): look at this; wrong constraint types!!
    session = interface.pymol.new_session()
    ligand_area:float=interface.pymol.general_cmd(session, [
        ('load', pdb_file),
        ('get_area', f"chain {ligand.parent.name} and resi {ligand.idx}")
    ])[-1]
    sasa_cutoff:int = int(max_sasa_ratio*ligand_area)

    options_file:str = create_docking_options(
                work_dir,
                pdb_file,
                xml_script,
                param_files,
                rng_seed,
                n_struct,
                ligand.parent.name,
                sasa_cutoff,
                clash_cutoff,
                contact_threshold,
                f"{ligand.idx}{ligand.parent.name}",
                cst_energy,
                cst_file,
                box_size,
                move_distance,
                transform_angle,
                transform_cycles,
                transform_repeats,
                transform_temperature,
                fr_repeats,
                grid_width
                )

    _LOGGER.info("Beginning RosettaLigand geometry sampling step...")
    
    start_dir: str = os.getcwd()
    
    os.chdir(work_dir)
    
    try:
        interface.rosetta.run_rosetta_scripts([f"@{Path(options_file).name}"])
    except:
        pass


    os.chdir(start_dir)

    scores_file:str = f"{work_dir}/docking_scores.sc"
    
    fs.check_file_exists(scores_file, exit_script = False )

    df: pd.DataFrame = interface.rosetta.parse_score_file(scores_file)

    df['description'] = df.apply(lambda row: f"{work_dir}/{row.description}.pdb", axis=1)

    energy_key:str=None        
    if use_qm:
        evaluate_geometry_qm_energy(df, structure, cluster_distance)
        energy_key='qm_energy'
    else:
        df['combined_energy'] = df.total_score + df.cst_filter
        energy_key='combined_energy'

    infile:str=df.sort_values(by=energy_key).description.to_list()[0]

    _parser = PDBParser()
    ref_stru = _parser.get_structure(infile)
    stru_oper.update_residues(structure, ref_stru)

    for cst in constraints:
        cst.change_topology(structure)

    

def create_docking_xml(work_dir:str) -> str:
    
    script_name:str=f"{work_dir}/docking_script.xml"

    content:List[str]="""<ROSETTASCRIPTS>
	<SCORINGGRIDS ligand_chain="%%ligand_chain%%" width="%%grid_width%%" name="grid">
		<ClassicGrid grid_name="classic" weight="1.0"/>
	</SCORINGGRIDS>
	<RESIDUE_SELECTORS>
		<Index name="ligand" resnums="%%ligand_idx%%"/>
		<CloseContact name="ligand_active_site" residue_selector="ligand" contact_threshold="%%contact_threshold%%"/>
		<Not name="not_ligand_active_site" selector="ligand_active_site"/>
	</RESIDUE_SELECTORS>
	<SCOREFXNS>
		<ScoreFunction name="ligand_soft_rep" weights="ligand_soft_rep">
			<Reweight scoretype="coordinate_constraint" weight="1.0"/>
			<Reweight scoretype="atom_pair_constraint" weight="1.0"/>
			<Reweight scoretype="angle_constraint" weight="1.0"/>
			<Reweight scoretype="dihedral_constraint" weight="1.0"/>
			<Reweight scoretype="chainbreak" weight="1.0"/>
		</ScoreFunction>
		<ScoreFunction name="hard_rep" weights="ligand"/>
	</SCOREFXNS>
	<LIGAND_AREAS>
	</LIGAND_AREAS>
	<INTERFACE_BUILDERS>
	</INTERFACE_BUILDERS>
	<MOVEMAP_BUILDERS>
	</MOVEMAP_BUILDERS>
	<TASKOPERATIONS>
	</TASKOPERATIONS>
	<SIMPLE_METRICS>
	</SIMPLE_METRICS>
    <CONSTRAINT_GENERATORS>
        <FileConstraintGenerator name="add_cst" filename="%%cst_file%%" />
    </CONSTRAINT_GENERATORS>
	<MOVERS>

		<Transform name="dock" chain="%%ligand_chain%%" box_size="%%box_size%%" move_distance="%%move_distance%%" angle="%%transform_angle%%" cycles="%%transform_cycles%%" repeats="%%transform_repeats%%" temperature="%%transform_temperature%%" grid_set="grid" use_constraints="true" cst_fa_file="%%cst_file%%"/>
        
        <ConstraintSetMover name="add_csts" add_constraints="true" cst_file="%%cst_file%%" />

        <ClearConstraintsMover name="rm_csts" />

        <FastRelax name="frelax" scorefxn="hard_rep" cst_file="%%cst_file%%" repeats="%%fr_repeats%%">
			<MoveMap name="full_enzyme" bb="true" chi="true" jump="true">
				<ResidueSelector selector="ligand_active_site"     bb="true" chi="true" bondangle="true" />
				<ResidueSelector selector="not_ligand_active_site" bb="false" chi="true" bondangle="false"/>
			</MoveMap>
		</FastRelax>
	</MOVERS>

	<FILTERS>
		<SimpleMetricFilter name="clash_filter" comparison_type="lt_or_eq" cutoff="%%clash_cutoff%%" composite_action="any">
			<PerResidueClashMetric name="clash" residue_selector="ligand" residue_selector2="ligand_active_site"/>
		</SimpleMetricFilter>
        <SimpleMetricFilter name="sasa_filter" comparison_type="lt_or_eq" cutoff="%%sasa_cutoff%%" composite_action="any">
            <SasaMetric name="sasa" residue_selector="ligand" />
		</SimpleMetricFilter>

        <ConstraintScore name="cst_filter" constraint_generators="add_cst" threshold="%%cst_cutoff%%" />
	</FILTERS>

	<PROTOCOLS>
		<Add mover_name="dock"      />
		<Add filter="clash_filter"  />
        <Add filter="sasa_filter"   />
		<Add mover_name="frelax"    />
        <Add filter="cst_filter"    />
        <Add mover_name="rm_csts"   />
	</PROTOCOLS>
	<OUTPUT>
	</OUTPUT>
</ROSETTASCRIPTS>""".splitlines()
    
    fs.write_lines(script_name, content)
    
    return script_name

def create_minimization_xml(work_dir:str) -> str:

    script_name:str=f"{work_dir}/minimization_script.xml"

    content:str="""<ROSETTASCRIPTS>
	<RESIDUE_SELECTORS>
		<Index name="ligand" resnums="%%ligand_idx%%"/>
		<CloseContact name="ligand_active_site" residue_selector="ligand" contact_threshold="%%contact_threshold%%"/>
		<Not name="not_ligand_active_site" selector="ligand_active_site"/>
	</RESIDUE_SELECTORS>
	<SCOREFXNS>
		<ScoreFunction name="hard_rep" weights="ligand"/>
	</SCOREFXNS>
	<LIGAND_AREAS>
	</LIGAND_AREAS>
	<INTERFACE_BUILDERS>
	</INTERFACE_BUILDERS>
	<MOVEMAP_BUILDERS>
	</MOVEMAP_BUILDERS>
	<TASKOPERATIONS>
	</TASKOPERATIONS>
	<SIMPLE_METRICS>
	</SIMPLE_METRICS>
    <CONSTRAINT_GENERATORS>
        <FileConstraintGenerator name="add_cst" filename="%%cst_file%%" />
    </CONSTRAINT_GENERATORS>
	<MOVERS>
        <FastRelax name="frelax" scorefxn="hard_rep" cst_file="%%cst_file%%" repeats="%%fr_repeats%%" ramp_down_constraints="%%ramp_constraints%%">
			<MoveMap name="full_enzyme" bb="true" chi="true" jump="true">
				<ResidueSelector selector="ligand_active_site"     bb="true" chi="true" bondangle="true" />
				<ResidueSelector selector="not_ligand_active_site" bb="false" chi="true" bondangle="false"/>
			</MoveMap>
		</FastRelax>
	</MOVERS>

	<FILTERS>
	</FILTERS>

	<PROTOCOLS>
		<Add mover_name="frelax" />
	</PROTOCOLS>
	<OUTPUT>
	</OUTPUT>
</ROSETTASCRIPTS>""".splitlines()

    fs.write_lines(script_name, content)
    
    return script_name

def create_minimization_options(work_dir:str,
                            pdb_file: str,
                            xml_file: str,
                            param_files: List[str],
                            rng_seed: int,
                            sele,
                            contact_threshold,
                            ramp_constraints,
                            cst_file,
                            fr_repeats
                            ) -> str:
                            
    """Makes the docking_options.txt file that the docking run will actually use. This function DOES NOT make any checks to the inputs.
    
    Args:
        pdb_file: The .pdb file (with constraints) to use. 
        xml_file: The validated RosettaScripts .xml file to be used for docking.
        param_files: The list() of reactant .params files.
        work_dir: The working directory 
        rng_seed: rng seed to be used during dcking.
        n_struct: Number of strutures to make as an int().

    Returns:
        Path to the docking_options.txt file with all the Rosetta options.
    """
    _LOGGER.info("Beginning creation of options file for Rosetta...")
    content: List[str] = [
        "-keep_input_protonation_state", 
        "-auto_setup_metals",
        "-run:constant_seed",
        f"-run:jran {int(rng_seed)}",
        "-in:file",
        f"    -s '{Path(pdb_file).name}'",
    ]

    for pf in param_files:
        content.append(f"    -extra_res_fa '{Path(pf).absolute()}'")

    stub_parent: str = os.path.expandvars(
        f"${config['rosetta.ROSETTA3']}/database/chemical/residue_type_sets/fa_standard/residue_types/protonation_states/")
    for stub in "GLU_P1.params GLU_P2.params LYS_D.params ASP_P1.params TYR_D.params HIS_P.params ASP_P2.params".split():
        content.append(f"    -extra_res_fa '{stub_parent}/{stub}'")

    content.extend([
        "-run:preserve_header",
        "-packing",
        "    -ex1",
        "    -ex2aro",
        "    -ex2 ",
        "    -no_optH false",
        "    -flip_HNQ true",
        "    -ignore_ligand_chi true",
        "-parser",
        f"   -protocol {Path(xml_file).absolute()}",
        f"""   -script_vars ligand_idx={sele} \\""",
        f"""                contact_threshold={contact_threshold} \\""",
        f"""                ramp_constraints={"true" if ramp_constraints else "false"} \\""",
        f"""                cst_file={cst_file} \\""",
        f"""                fr_repeats={fr_repeats}""",
        "-out",
        f"   -file:scorefile 'minimization_scores.sc'",
        "   -level 200",
        "   -overwrite",
        "   -path",
        f"       -all '{work_dir}'",
    ])

    fname = Path(work_dir) / "minimization_options.txt"
    score_file: str = f"{work_dir}/minimization_scores.sc"

    _LOGGER.info(f"\toptions file: {fname}")
    _LOGGER.info(f"\tscore file: {score_file}")
    _LOGGER.info(f"\tenzyme-reactant complexes directory: {work_dir}/complexes")
    
    fs.safe_rm(fname)
    fs.safe_rm(score_file)
    
    _LOGGER.info(f"Wrote the below settings to {fname}:")
    for ll in content:
        _LOGGER.info(f"\t{ll}")
    fs.write_lines(fname, content)

    option_file = fname.absolute()

    return str(fname)



def create_docking_options(work_dir:str,
                            pdb_file: str,
                            xml_file: str,
                            param_files: List[str],
                            rng_seed: int,
                            n_struct: int,
                            ligand_chain:str,
                            sasa_cutoff:int,
                            clash_cutoff:int,
                            contact_threshold:float,
                            ligand_idx:str,
                            cst_cutoff:int,
                            cst_file:str,
                            box_size:float,
                            move_distance:float,
                            transform_angle:int,
                            transform_cycles:int,
                            transform_repeats:int,
                            transform_temperature:int,
                            fr_repeats:int,
                            grid_width:float) -> str:
                            
    """Makes the docking_options.txt file that the docking run will actually use. This function DOES NOT make any checks to the inputs.
    
    Args:
        pdb_file: The .pdb file (with constraints) to use. 
        xml_file: The validated RosettaScripts .xml file to be used for docking.
        param_files: The list() of reactant .params files.
        work_dir: The working directory 
        rng_seed: rng seed to be used during dcking.
        n_struct: Number of strutures to make as an int().

    Returns:
        Path to the docking_options.txt file with all the Rosetta options.
    """
    _LOGGER.info("Beginning creation of options file for Rosetta...")
    content: List[str] = [
        "-keep_input_protonation_state", 
        "-auto_setup_metals",
        "-run:constant_seed",
        f"-run:jran {int(rng_seed)}",
        "-in:file",
        f"    -s '{Path(pdb_file).name}'",
    ]

    for pf in param_files:
        content.append(f"    -extra_res_fa '{Path(pf).absolute()}'")

    stub_parent: str = os.path.expandvars(
        f"${config['rosetta.ROSETTA3']}/database/chemical/residue_type_sets/fa_standard/residue_types/protonation_states/")
    for stub in "GLU_P1.params GLU_P2.params LYS_D.params ASP_P1.params TYR_D.params HIS_P.params ASP_P2.params".split():
        content.append(f"    -extra_res_fa '{stub_parent}/{stub}'")

    content.extend([
        "-run:preserve_header",
        "-packing",
        "    -ex1",
        "    -ex2aro",
        "    -ex2 ",
        "    -no_optH false",
        "    -flip_HNQ true",
        "    -ignore_ligand_chi true",
        "-parser",
        f"   -protocol {Path(xml_file).absolute()}",
        f"""   -script_vars ligand_chain={ligand_chain} \\""",
        f"""                {sasa_cutoff=} \\""",
        f"""                {clash_cutoff=} \\""",
        f"""                contact_threshold={contact_threshold:.1f} \\""",
        f"""                ligand_idx={ligand_idx} \\""",
        f"""                {cst_cutoff=} \\""",
        f"""                cst_file={cst_file} \\""",
        f"""                box_size={box_size:.1f} \\""",
        f"""                move_distance={move_distance:.1f} \\""",
        f"""                {transform_angle=} \\""",
        f"""                {transform_cycles=} \\""",
        f"""                {transform_repeats=} \\""",
        f"""                {transform_temperature=} \\""",
        f"""                {fr_repeats=} \\""",
        f"""                grid_width={grid_width:.1f}""", 
        "-out",
        f"   -file:scorefile 'docking_scores.sc'",
        "   -level 200",
        f"   -nstruct {n_struct}",
        "   -overwrite",
        "   -path",
        f"       -all '{work_dir}'",
    ])

    qsar_grid: str = str(Path(f"{work_dir}/qsar_grids/").absolute())
    content.append(f"-qsar:grid_dir {qsar_grid}")

    fname = Path(work_dir) / "docking_options.txt"
    score_file: str = f"{work_dir}/docking_scores.sc"

    _LOGGER.info(f"\toptions file: {fname}")
    _LOGGER.info(f"\tscore file: {score_file}")
    _LOGGER.info(f"\tenzyme-reactant complexes directory: {work_dir}/")
    _LOGGER.info(f"\tqsar_gird directory: {qsar_grid}")
    
    fs.safe_rm(fname)
    fs.safe_rm(score_file)
    fs.safe_rmdir(qsar_grid)
    fs.safe_mkdir(qsar_grid)
    
    _LOGGER.info(f"Wrote the below settings to {fname}:")
    for ll in content:
        _LOGGER.info(f"\t{ll}")
    fs.write_lines(fname, content)

    option_file = fname.absolute()

    return str(fname)


def parameterize_structure(structure:Structure, work_dir:str) -> List[str]:
    
    result:List[str] = list()
    for residue in structure.residues:
        if residue.is_ligand():
            result.append(generate_ligand_params(residue, work_dir))

    return result

def generate_ligand_params(ligand:Ligand, work_dir:str) -> List[str]:
    """Given the input Structure(), parameterize everything needed to use the Ligand()'s in Rosetta.
    
    Args:
        stru: The Structure() to parameterize.
        work_dir: Where temporary files will be saved.

    Returns:
        A List[str] with ligand .params files.         
    """
    #TODO(CJ): probably move this to the RosettaInterface
    _LOGGER.info("Beginning preparation of each reactant...")
  
    param_file:str = f"{work_dir}/{ligand.name}.params"
    if Path(param_file).exists():
        return param_file


    parser = Mol2Parser()
    if ligand.net_charge is None:
        ligand.net_charge = interface.bcl.calculate_formal_charge( ligand )


    param_file:str = interface.rosetta.parameterize_ligand(ligand, charge=ligand.net_charge, work_dir=work_dir)
    
    _LOGGER.info(f"Information for reactant {ligand.name}:")
    _LOGGER.info(f"\tparam file: {param_file}")
    _LOGGER.info(f"\tcharge: {ligand.net_charge}")
    _LOGGER.info(f"\tnum conformers: {ligand.n_conformers()}")
    
    _LOGGER.info("Finished reactant preparation!")
    return param_file


def get_active_site_sele(structure: Structure, distance_cutoff:float) -> str:
    """Creates a pymol-compatible sele for the active site of the supplied Structure(). Basic structure
    is to select all residues within the specified cutoff of the Ligand()'s in the structure. Metal ions up to
    2*distance_cutoff from the Ligand()'s are also selected.
    
    Args:
        structure: The Structure() to be used as a template for the active site.
        distance_cutoff: The cutoff in Angstroms for a Residue() to be included in the active site.  
        
    Returns:
        Selection string in pymol format which defines the enzyme active site.        
    """

    _LOGGER.info("Analyzing enzyme active site...")
    #for res in structure.residues:
    ligand_residue_keys = set()
    for res in structure.residues:
        if res.is_ligand():
            ligand_residue_keys.add(res.key())

    parser = PDBParser()
    start_pdb:str=f"{config['system.SCRATCH_DIR']}/active_site_selection.pdb"
    parser.save_structure(start_pdb, structure)
   
    ligand_sele:str = f"byres ( {'or '.join(map(lambda lrk: f'(all within {distance_cutoff} of chain {lrk[0]} and resi {lrk[1]})', ligand_residue_keys))}  )"
    ligand_sele:str = " or ".join(map(
        lambda lrk: f"((byres all within {distance_cutoff:.2f} of (chain {lrk[0]} and resi {lrk[1]})) or (metals within  {2*distance_cutoff:.2f} of( chain {lrk[0]} and resi {lrk[1]})))",
        ligand_residue_keys
    ))

    session = interface.pymol.new_session()
    df = interface.pymol.collect(session, start_pdb, "chain resi".split(), sele=ligand_sele)

    fs.safe_rm(start_pdb)

    result = set()
    for i, row in df.iterrows():
        result.add((row['chain'], row['resi']))
    
    _LOGGER.info(f"Found {len(result)} residues within {distance_cutoff} angstroms of reactants!")
    
    return " or ".join(map(
        lambda rr: f"( chain {rr[0]} and resi {rr[1]})",
        result
    ))

def qm_minimization(structure:Structure,
                constraints:List[StructureConstraint],
                cluster_distance:float,
                freeze_ligands:bool,
                work_dir:str) -> None:
    """Performs QM minimization of the enzyme active site using xtb. Assumes that backbone atoms of the Residue()'s should be 
    frozen. Is capable of converting supplied constraints to xtb format. Updates coordinates in place. 
    
    Args:
        structure: The Structure() object to optimize.
        constraints: The List[StructureConstraint] which define the enzyme active site geometry.
        cluster_distance: The cutoff in angstroms for inclusion of a Residue() in the active site.
        work_dir: The temporary directory to do all work.

    Returns:
        Nothing.
    """

    as_sele = get_active_site_sele(structure, cluster_distance)
   
    to_freeze:List[Atom] = list()

    for res in structure.residues:
        if res.is_canonical():
            continue

        for atom in res.atoms:
            atom.charge = 0.0

            if freeze_ligands and atom.element != 'H':
                to_freeze.append(atom)
    
    translate_structure(structure, start_naming='rosetta')
    es = qm_optimize(structure,
            engine="xtb",
            constraints=constraints + [CartesianFreeze(structure.backbone_atoms() + to_freeze )],
            regions=[as_sele],
            region_methods=[chem.QMLevelOfTheory(basis_set='',method='GFN2', solvent='water', solv_method='ALPB')],
            parallel_method=None)[0]

    translate_structure(structure, end_naming='rosetta')

def evaluate_geometry_qm_energy(df: pd.DataFrame, structure: Structure, cluster_cutoff: float) -> None:
    """Aids in ranking and selection of candidate geometries through a semi-empirical QM single point energy
    calculation with xtb. Creates a capped active site by using a specified cluster_cutoff parameter to specify
    the enzyme's active site.

    Args:
        df: The geometry DataFrame containing all information  
        structure: The reference Structure() in use.
        cluster_cutoff: The cutoff in Angstroms for a Residue() to be included in the QM region. 

    Returns:
        Nothing.        
    """
#    if not df.selected.sum():
#        _LOGGER.error("No geometries are still selected!")
#        raise TypeError()

    _LOGGER.info(f"Beginning qm energy evaluation.")
    as_sele:str = get_active_site_sele(structure, cluster_cutoff)
    qm_energy = []

    _parser = PDBParser()
    for i, row in df.iterrows():


        energy:float = None
        _df_stru = _parser.get_structure( row.description )
        translate_structure(_df_stru, start_naming='rosetta')
        for res in structure.residues:
            if res.is_canonical():
                continue
            _df_stru.get(res.key_str).net_charge = res.net_charge
            _df_stru.get(res.key_str).multiplicity = res.multiplicity
            for atom in _df_stru.get(res.key_str):
                atom.charge = 0.0
        es = single_point(
            _df_stru,
            engine='xtb',
            region_methods=[chem.QMLevelOfTheory(basis_set='',method='GFN2', solvent='water', solv_method='ALPB')],
            parallel_method=None,
            regions=[as_sele])
        qm_energy.append(es[0].energy_0)

    df['qm_energy'] = qm_energy
    _LOGGER.info("Finished qm energy evaluation!")
