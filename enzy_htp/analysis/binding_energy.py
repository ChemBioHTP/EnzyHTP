"""
"""
from enzy_htp.core import file_system as fs
from pathlib import Path
from enzy_htp import interface, config

from enzy_htp.structure import PDBParser, translate_structure 

from typing import Dict, Callable
#from enzy_htp._interface.rosetta_interface import (
#        RosettaOptions,
#        RosettaScriptsElement,
#        RosettaScriptsProtocol)

from enzy_htp.structure import Structure, Atom, Ligand

def binding_energy(
    stru: Structure,
    ligand: Ligand,
    engine:str,
    rosetta_sfxn:"RosettaScriptsElement"=None,
    work_dir:str=None
    ) -> float:

    
    if engine == 'rosetta':
        return binding_energy_rosetta(stru, ligand, rosetta_sfxn, work_dir)





def binding_energy_rosetta( stru: Structure, ligand: Ligand, rosetta_sfxn:"RosettaScriptsElement", work_dir) -> float:


    if work_dir is None:
        work_dir = "./"

    fs.safe_mkdir( work_dir )

    if rosetta_sfxn is None:

        rosetta_sfxn = RosettaScriptsElement("ScoreFunction", name="ligand_soft_rep", weights="ligand_soft_rep",
            children=[
                RosettaScriptsElement("Reweight", scoretype="fa_elec", weight="0.42"),
                RosettaScriptsElement("Reweight", scoretype="hbond_bb_sc", weight="1.3"),
                RosettaScriptsElement("Reweight", scoretype="hbond_sc", weight="1.3"),
                RosettaScriptsElement("Reweight", scoretype="rama", weight="0.2"),
        ])

    translate_structure(stru, end_naming='rosetta')

    chain_names = set()
    for res in stru.residues:
        chain_names.add(res.parent.name)
        if res.key() == ligand.key():
            break
    else:
        assert False

    jump:int = len(chain_names) - 1 

    fname:str=f"{work_dir}/binding_energy_script.xml"
    
    #TODO(CJ): put in the part about where the ligands will get parameterized
    
    protocol = RosettaScriptsProtocol()

    protocol.add_scorefunction( rosetta_sfxn )
    protocol.add_filter( RosettaScriptsElement("Ddg", name="binding_dg", jump=f"{jump}", threshold="100", scorefxn=f"{rosetta_sfxn.attrib['name']}"))
    protocol.add_protocol( RosettaScriptsElement("Add", filter_name="binding_dg") )

    parser = PDBParser()    
    pdb_file = f"{work_dir}/rosetta_start.pdb"
    parser.save_structure(pdb_file, stru)

    opts = RosettaOptions() 
    opts['in:file:s'] = str(Path(pdb_file).absolute())
    opts['overwrite'] = True
    opts['out:file:scorefile'] = str(Path(f"{work_dir}/score.sc").absolute())

    fs.safe_rm(opts['out:file:scorefile'])

    rs_engine = interface.rosetta.build_rosetta_scripts_engine( protocol, opts, work_dir )



    updated_structure = rs_engine.run(stru)[0]
    return updated_structure.data['binding_dg']

