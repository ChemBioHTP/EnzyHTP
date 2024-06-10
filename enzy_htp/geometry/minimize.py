

from pathlib import Path
from enzy_htp import interface
from enzy_htp.core import file_system as fs
from typing import List, Dict
from enzy_htp.structure import (
    Structure,
    StructureEnsemble,
    StructureConstraint,
    PDBParser,
    Mol2Parser,
    structure_translator,
)

from enzy_htp.structure.structure_operation import    update_residues

from enzy_htp._interface.rosetta_interface import (
    RosettaOptions,
    RosettaScriptsElement,
    RosettaScriptsProtocol
)



def minimize( 
    stru:Structure,
    method='rosetta',
    work_dir:str="./minimize",
    constraints:List[StructureConstraint]=None,
    movemap:List[Dict]=None,
    n_iter:int=5,
    n_struct:int=1,
    **kwargs
):
    
    fs.safe_mkdir( work_dir )

    if method == 'rosetta':
        return minimize_rosetta(stru, n_iter, n_struct, movemap, work_dir, **kwargs)


def minimize_rosetta(stru, n_iter, n_struct, movemap, work_dir, **kwargs):
    
    params = list()

    structure_translator.translate_structure(stru, end_naming='rosetta')

    for res in stru.residues:
        if res.is_ligand():
            params.append(Path(interface.rosetta.parameterize_ligand(res, work_dir=work_dir)).absolute())

    print(params)
    options = RosettaOptions()
    options['extra_res_fa'] = params
    options['nstruct'] = n_struct
    options['overwrite'] = 'true'
    options['packing:no_optH'] = 'false'

    protocol = RosettaScriptsProtocol()
    protocol.add_scorefunction(RosettaScriptsElement("ScoreFunction", name="r15", weights="ref2015"))
    protocol.add_protocol(RosettaScriptsElement("Add", mover_name="frelax"))
    mm_children = list()

    temp_pdb = f"{work_dir}/min_temp.pdb"
    parser = PDBParser()
    parser.save_structure(temp_pdb, stru)
    
    for midx,mm in enumerate(movemap):
        session = interface.pymol.new_session() 
        df = interface.pymol.collect(session, temp_pdb, sele=mm['sele'], variables='chain resi'.split())
        rs=",".join(map(lambda pr: pr[-1]+pr[0],
            set(list(zip(df.chain,df.resi)))
        ))
        
        rs_name=f"res_selector_{midx}"
        protocol.add_residue_selector(RosettaScriptsElement("Index", name=rs_name, resnums=rs))
        
        mm_children.append(RosettaScriptsElement("ResidueSelector", selector=rs_name, bb=mm.get('bb', 'false'), chi=mm.get('chi', 'false'), bondangle=mm.get('bondangle', 'false')))

    fs.safe_rm(temp_pdb)

    protocol.add_mover(RosettaScriptsElement("FastRelax", name="frelax", scorefxn="r15", repeats=str(n_iter),
        children=[RosettaScriptsElement("MoveMap", name="full_structure", bb="true", chi="true", jump="true", children=mm_children)]
    ))
    score_file = interface.rosetta.run_rscripts(stru, protocol, options, work_dir)

    
    df = interface.rosetta.parse_score_file(score_file)
    new_fname = f"{work_dir}/{df.iloc[0].description}.pdb"
    stru2 = parser.get_structure(new_fname)
    update_residues(stru, stru2)
#    print(set([aa.name for aa in stru2.residues[2].atoms]) == set([aa.name for aa in stru.residues[2].atoms]))
#    for aa in stru.atoms:
#        aa.coord = stru2.get(aa.key).coord
#
#    structure_translator.translate_structure(stru, start_naming='rosetta')
