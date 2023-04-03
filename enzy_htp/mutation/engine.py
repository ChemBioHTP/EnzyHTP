"""Top level engine for mutating the structure in a .pdb file. Users should only interface with the mutate_pdb()
method which gives users the ability to mutate a .pdb file and specify all relevant aspects of this process.
Mutation is carried out by an underlying engine and the supported engines currently include:
    + Amber/tleap
    + PyMOL

Note that the current implementation will mutate the .pdb file and keep the residue indicies and chain names 
consistent.
Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-15
"""

#TODO(CJ): This doesn't currently work for tleap. Get it working

import hashlib
from pathlib import Path
from copy import deepcopy
from string import ascii_uppercase
from collections import defaultdict
from typing import List, Dict, Union, Any

import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb

from enzy_htp import interface
import enzy_htp.chemical as chem
import enzy_htp.structure as struct
import enzy_htp.preparation as prep
from enzy_htp.core import file_system as fs
from enzy_htp.core import _LOGGER, UnsupportedMethod

from .mutation_restrictions import MutationRestrictions, restriction_object
from .mutation import generate_all_mutations, Mutation


def mutate_pdb(
    pdb: str,
    n_mutations: int = 1,
    mutations: List[Mutation] = list(),
    restrictions: MutationRestrictions = None,
    engine: str = "tleap",
    out_dir: str = None,
    random_state: int = 100,
) -> str:
    """Top level method that users should call to mutated a given .pdb file. Function gives the ability to
    specify the number of mutations, specific mutations to employ, the mutation engine and the random state
    of the random number generator. Further information and additional parameters is provided in the Args
    section below. Information about the specific mutations applied to the structure are encoded into the
    new .pdb file. The new .pdb file is returned from the main function and has the format:
            <base_name>_<mut1>_<mut2>...<mutN>.pdb
    Where <mutN> has the format:
        <original one-letter name><residue number><mutated one-letter name>
    A commandline example:
    >>> mutated_pdb:str = mutate_pdb("enzy.pdb")
    >>> mutated_pdb
    "enzy_A10G.pdb"

    Args:
        pdb: The name of the base .pdb file. Only required value.
        n_mutations: The number of desired mutations as an int(). Default value is 1.
        mutations: A list() of Mutation namedtuple()'s. Default value is None.
        restrictions: A MutationRestrictions() object to restrict generated mutations. Default is None.
        engine: The name of the mutation engine as a str().
        out_dir: The directory to save the mutated .pdb to. Default value is None, in which case
            the file is saved to the current directory. When it doesn't exist, out_dir is created.
        random_state: The int() seed for the random number generator. Default value is 100.
    Raises:
        enzy_htp.core.exception.UnsupportedMethod if the supplied engine is not supported.
    Returns:
        The name of the .pdb file with the mutated structure.
    """

    def unique_by_pos(muts) -> List[Mutation]:
        holder = defaultdict(list)
        for mm in muts:
            holder[(mm.chain_id, mm.res_num)].append(deepcopy(mm))
        result = list()
        for mlist in holder.values():
            result.append(mlist.pop())
        return result

    mutations = unique_by_pos(mutations)
    needed: int = n_mutations - len(mutations)
    if needed > 0:
        mutations.extend(get_mutations(pdb, needed, mutations, random_state,
                                       restrictions))
    else:
        _LOGGER.warning(
            f"{len(mutations)} unique mutations were provided but only {n_mutations} were requested."
        )
        _LOGGER.warning(
            f"The first {n_mutations} will be used. To change this behavior, increase n_mutations."
        )
        mutations = sorted(mutations)[0:n_mutations]

    if out_dir is not None:
        fs.safe_mkdir(out_dir)

    outfile: str = mutated_name(pdb, out_dir, mutations)

    IMPLEMENTATION: Dict = {
        "tleap": _mutate_tleap, 
        "pymol": _mutate_pymol,
        "rosetta": _mutate_rosetta
        }

    if engine not in IMPLEMENTATION:
        raise UnsupportedMethod(
            f"'{engine}' is not a supported mutation engine. Allowed values are {', '.join(list(IMPLEMENTATION.keys()))}"
        )

    IMPLEMENTATION[engine](pdb, outfile, mutations)

    return outfile


def get_mutations(
    pdb: str,
    needed: int,
    mutations: List[Mutation],
    random_state: int,
    restrictions: Union[MutationRestrictions, None],
) -> List[Mutation]:
    """Method that curates Mutation() namedtuples given a number of restrictions including
    a target number of mutations to get, the existing mutations and the MutationRestrictions() object.
    Note that this method SHOULD NOT be called directly by the user.
    Args:
        pdb: A str() with the file corresponding to the original structure.
        needed: The number of desired Mutation() objects as an int().
        mutations: A list() of existing Mutation() namedtuples.
        random_state: The int() rng seed for the rng engine.
        restrictions: The MutationRestrictions() object. Can be None.
    Raises:
        A base Exception() if not enough Mutation namedtuple()'s can be generated.
    Returns:
        A list() of Mutation() namedtuples compatible with the supplied enzyme system and constraints.
    """

    def random_list_elem(ll: list) -> Any:
        """Helper method that randomly chooses an element from a list. numpy.random.choice() doesn't 
        like to take elements from list()'s of tuples so this is the work around.
        """
        return ll[np.random.randint(len(ll))]

    if restrictions is None:
        restrictions = restriction_object(pdb)

    structure: struct.Structure = struct.PDBParser.get_structure(pdb)
    mut_dict: Dict[Tuple[Str, int], List[Mutation]] = generate_all_mutations(structure)

    for mut in mutations:
        restrictions.lock_residue((mut.chain_id, mut.res_num))

    mut_dict = restrictions.apply(mut_dict)

    result: List[Mutation] = list()

    np.random.seed(random_state)

    while len(result) < needed:

        if not len(mut_dict):
            raise Exception(
                "Could not generate enough Mutation namedtuple()'s. This is likely due to too many restrictions being present."
            )

        md_keys = list(mut_dict.keys())
        chosen = random_list_elem(md_keys)
        if len(mut_dict[chosen]):
            result.append(random_list_elem(sorted(mut_dict[chosen])))
        del mut_dict[chosen]

    return result


def mutated_name(pdb: str, outdir: str, mutations: List[Mutation]) -> str:
    """Helper method that generates a PDB file name which encodes information about the
    mutations which will be applied to the structure in the file. In the case that the original residue
    is listed as 'X' for unknown, the identity will be determined from the input .pdb file.
    Args:
        pdb: Filename for the .pdb as a str().
        outdir: The directory to output the mutated .pdb. If None, uses the same base directory as the pdb parameter.
        mutations: A list() of the Mutation() namedtuples to be applied to the structure.
    Returns:
        The filename that the mutated structure should be saved to.
    """
    pdb_df: pd.DataFrame = PandasPdb().read_pdb(pdb).df["ATOM"]

    def deduce_original(mm: Mutation, df: pd.DataFrame) -> str:
        for i, row in df.iterrows():
            if row.chain_id == mm.chain_id and row.residue_number == mm.res_num:
                return chem.convert_to_one_letter(row.residue_name)
        else:
            # TODO(CJ): improve this
            raise TypeError()

    pdb_path = Path(pdb)
    if outdir is None:
        outdir = str(pdb_path.parent)

    name_stem: str = pdb_path.stem
    mutations = sorted(mutations, key=lambda m: (m.chain_id, m.res_num))
    for mut in mutations:
        if mut.orig == "X":
            orig_name: str = deduce_original(mut, pdb_df)
        else:
            orig_name: str = mut.orig
        name_stem += f"_{orig_name}{mut.res_num}{mut.target}"
    return f"{outdir}/{name_stem}.pdb"


def _mutate_pymol(pdb: str, outfile: str, mutations: List[Mutation]) -> None:
    """Underlying implementation of mutation with pymol. Serves as implementation only, SHOULD NOT
    BE CALLED BY USERS DIRECTLY. Follows generalized function signature taking the name of the .pdb,
    the outfile to save the mutated version to and a list() of mutations. Function assumes that the 
    supplied mutations are valid. Procedure utilizes the cmd.wizard() api exposed by the pymol python 
    package.
    
    Args:
        pdb: The name of the original .pdb file as a str().
        outfile: Name of the file to save the mutated structure to.
        mutations: A list() of Mutation namedtuple()'s to apply.
    Returns:
        Nothing.

    """
    #TODO(CJ): finish this
    pass


def _mutate_rosetta(pdb:str, outfile:str, mutations: List[Mutation]) -> None:
    """Underlying implementation of mutation with Rosetta. Serves as an implementation only,
    SHOULD NOT BE CALLED BY USERS DIRECTLY. Follows generalized function signature taking the 
    name of the .pdb, the oufile to save the mutated version to and a list() of mutations. 
    Function assumes that the supplied mutations are valid. Procedure is to perform the mutation
    using RosettaScripts and the MutateResidue mover.
    Args:
        pdb: The name of the original .pdb file as a str().
        outfile: Name of the file to save the mutated structure to.
        mutations: A list() of Mutation namedtuple()'s to apply.
    Returns:
        Nothing.

    """

    working_dir = Path(pdb).parent

    contents:List[str] = ["""<ROSETTASCRIPTS>
    <SCOREFXNS>
    </SCOREFXNS>
    <RESIDUE_SELECTORS>
    </RESIDUE_SELECTORS>
    <TASKOPERATIONS>
    </TASKOPERATIONS>
    <SIMPLE_METRICS>
    </SIMPLE_METRICS>
    <FILTERS>
    </FILTERS>
    <MOVERS>
"""]
    mut_names:List[str] = list()
    
    #TODO(CJ): have to solve the rosetta indexing problem
    for midx, mut in enumerate(mutations):
        index:int = resolve_index(pdb_lines, mut)
        mover_name:str=f"mr{midx+1}"
        contents.append(f"        <MutateResidue name=\"{mover_name}\" target=\"{TODO}\" new_res=\"{TODO}\" />")
        mut_names.append( mover_name )

    contents.append("""    </MOVERS>
    <PROTOCOLS>
""")
    
    for mn in mut_names:
        contents.append(f"        <Add mover_name=\"{mn}\" />")
    contents.append("""    </PROTOCOLS>
    <OUTPUT />
<ROSETTASCRIPTS>
""")
    
    script_file:Path = working_dir / "enzy_htp_mut.xml"
    fs.write_lines(script_files, contents)
   

    options:List[str] = list()
    interface.rosetta.run_rosetta_scripts(
        pdb,
        script_file,
        options
    )

    #TODO(CJ): have to check that the file exists in the working directory

    #TODO(CJ): copy the file to the output file


def _mutate_tleap(pdb: str, outfile: str, mutations: List[Mutation]) -> None:
    """Underlying implementation of mutation with tleap. Serves as implementation only, SHOULD NOT
    BE CALLED BY USERS DIRECTLY. Follows generalized function signature taking the name of the .pdb,
    the oufile to save the mutated version to and a list() of mutations. Function assumes that the
    supplied mutations are valid. Procedure is to replace backbone atoms and their names in the .pdb
    file and then apply the changes with Amber's tleap tool.
    Args:
        pdb: The name of the original .pdb file as a str().
        outfile: Name of the file to save the mutated structure to.
        mutations: A list() of Mutation namedtuple()'s to apply.
    Returns:
        Nothing.
    """

    def get_backup(fname: str):
        lines = prep.read_pdb_lines(fname)
        lines = list(filter(lambda ll: ll.is_ATOM() or ll.is_HETATM(), lines))
        holder = defaultdict(list)
        for ll in lines:
            holder[f"{ll.chain_id}.{ll.resi_name}.{ll.resi_id}"].append(ll)
        to_remove = list()
        for key in holder.keys():
            if key.split('.')[1] in chem.THREE_LETTER_AA_MAPPER:
                to_remove.append(key)
        for tr in to_remove:
            del holder[tr]
        return holder

    backup = get_backup(pdb)
    chain_count: int = 1
    chain_names: List[str] = list(ascii_uppercase[::-1])
    curr_chain: str = chain_names.pop()
    pdb_lines: List[prep.PDBLine] = prep.read_pdb_lines(pdb)
    mask: List[bool] = [True] * len(pdb_lines)
    for pidx, pdb_l in enumerate(pdb_lines):
        if pdb_l.is_TER():
            curr_chain = chain_names.pop()
            continue
        # only match in the dataline and keep all non data lines
        if not pdb_l.is_ATOM():
            continue
        for mf in mutations:
            # Test for every Flag for every lines

            if curr_chain == mf.chain_id and pdb_l.resi_id == mf.res_num:
                # fix for mutations of Gly & Pro
                old_atoms: List[str] = {
                    "G": ["N", "H", "CA", "C", "O"],
                    "P": ["N", "CA", "HA", "CB", "C", "O"],
                }.get(mf.target, ["N", "H", "CA", "HA", "CB", "C", "O"])

                line: str = pdb_l.line
                for oa in old_atoms:
                    if oa == pdb_l.atom_name:
                        pdb_l.line = f"{line[:17]}{chem.convert_to_three_letter(mf.target)}{line[20:]}"
                        break
                else:
                    mask[pidx] = False

    fs.write_lines(outfile, np.array(list(map(lambda pl: pl.line, pdb_lines)))[mask])

    interface.amber.mutate(outfile)
    if backup:
        structure: struct.Structure = struct.structure_from_pdb(outfile)
        for rkey in structure.residue_keys:
            if rkey not in backup:
                continue
            tmp_file = f"/tmp/ligand.{rkey}.tmp.pdb"
            fs.write_lines(tmp_file,
                           list(map(lambda ll: ll.line, backup[rkey])) + ['END'])
            lig: struct.Ligand = struct.ligand_from_pdb(tmp_file)
            (cname, rname, rnum) = rkey.split('.')
            rnum = int(rnum)
            structure.chain_mapper[cname]._residues[rnum - 1] = lig
            fs.safe_rm(tmp_file)

        structure.to_pdb(outfile)
