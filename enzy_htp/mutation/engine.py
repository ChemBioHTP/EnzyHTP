"""Top level engine for mutating the structure in a .pdb file.


Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-06-15
"""
import numpy as np
from pathlib import Path
from string import ascii_uppercase
from typing import List, Dict, Union


import enzy_htp.chemical as chem
from enzy_htp.core import _LOGGER
from enzy_htp.core import file_system as fs
import enzy_htp.structure as struct
import enzy_htp.molecular_mechanics as mm
import enzy_htp.preparation as prep
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
    """TODO(CJ)"""
    # quickly remove exact duplicates
    mutations = list(set(mutations))
    needed: int = n_mutations - len(mutations)
    if needed > 0:
        mutations.extend(
            get_mutations(pdb, needed, mutations, random_state, restrictions)
        )

    outfile: str = mutated_name(pdb, out_dir, mutations)

    IMPLEMENTATION: Dict = {"tleap": _mutate_tleap}

    if engine not in IMPLEMENTATION:
        raise TypeError()

    IMPLEMENTATION[engine](pdb, outfile, mutations)


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

    Returns:
        A list() of Mutation() namedtuples compatible with the supplied enzyme system and constraints.
    """
    if restrictions is None:
        restrictions = restriction_object(pdb)

    structure: struct.Structure = struct.structure_from_pdb(pdb)
    mut_dict: Dict[Tuple[Str, int], List[Mutation]] = generate_all_mutations(structure)

    for mut in mutations:
        restrictions.lock_residue((mut.chain_id, mut.res_num))

    mut_dict = restrictions.apply(mut_dict)

    result: List[Mutation] = list()

    np.random.seed(random_state)

    while len(result) < needed:

        if not len(mut_dict):
            raise TypeError()

        md_keys = list(mut_dict.keys())
        chosen = np.random.choice(md_keys)
        if len(mut_dict[chosen]):
            result.append(np.random.choice(mut_dic[chosen]))
        del mut_dict[chosen]

    return result


def mutated_name(pdb: str, outdir: str, mutations: List[str]) -> str:
    """Helper method that generates a PDB file name which encodes information about the
    mutations which will be applied to the structure in the file.

    Args:
        pdb: Filename for the .pdb as a str().
        outdir: The directory to output the mutated .pdb. If None, uses the same base directory as the pdb parameter.
        mutations: A list() of the Mutation() namedtuples to be applied to the structure.

    Returns:
        The filename that the mutated structure should be saved to.
    """

    pdb_path = Path(pdb)
    if outdir is None:
        outdir = str(pdb_path.parent)

    name_stem: str = pdb_path.stem
    mutations = sorted(mutations, key=lambda m: (m.chain_id, m.res_num))
    for mut in mutations:
        name_stem += f"_{mut.orig}{mut.res_num}{mut.target}"
    return f"{outdir}/{name_stem}.pdb"


def _mutate_tleap(pdb: str, outfile: str, mutations: List[Mutation]) -> None:
    """Underlying implementation of mutation with tleap. Serves as impelementation only, SHOULD NOT
    BE CALLED BY USERS DIRECTLY. Follows generalized function signature taking the name of the .pdb,
    the oufile to save the mutated version to and a list() of mutations. Function assumes that the 
    supplied mutations are valid.

    Args:
        pdb: The name of the original .pdb file as a str().
        outfile: Name of the file to save the mutated structure to.
        mutations: A list() of Mutation namedtuple()'s to apply.  
    
    Returns:
        Nothing.
    """
    chain_count:int = 1
    chain_names:List[str] = list(ascii_uppercase[::-1])
    curr_chain:str = chain_names.pop()
    pdb_lines:List[prep.PDBLine] = prep.read_pdb_lines(pdb)
    mask:List[bool] = [True] * len(pdb_lines)
    for pdb_l in pdb_lines:
        if pdb_l.is_TER():
            curr_chain = chain_names.pop()
        match = 0
        # only match in the dataline and keep all non data lines
        if not pdb_l.is_ATOM():
            continue
        for mf in mutations:
            # Test for every Flag for every lines

            if (
                curr_chain == mf.chain_id
                and pdb_l.resi_id == mf.res_num
            ):
                # do not write old line if match a MutaFlag
                match = 1
                # fix for mutations of Gly & Pro
                old_atoms:List[str] = {
                    "G": ["N", "H", "CA", "C", "O"],
                    "P": ["N", "CA", "HA", "CB", "C", "O"],
                }.get(mf.target, ["N", "H", "CA", "HA", "CB", "C", "O"])

                line:str = pdb_l.line
                for oa in old_atoms:
                    if oa == pdb_l.atom_name:
                        pdb_l.line = f"{line[:17]}{chem.convert_to_three_letter(mf.target)}{line[20:]}"

    fs.write_lines(outfile, list(map(lambda pl: pl.line, pdb_lines)))
    ai = mm.AmberInterface()
    ai.mutate(outfile)
