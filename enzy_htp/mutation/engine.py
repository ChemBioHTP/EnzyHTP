"""Top level engine for mutating the structure in a .pdb file. Users should only interface with the mutate_pdb()
method which gives users the ability to mutate a .pdb file and specify all relevant aspects of this process.
Mutation is carried out by an underlying engine and the supported engines currently include:
    + Amber/tleap

Note that the current implementation will mutate the .pdb file and keep the residue indicies and chain names 
consistent.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>

Date: 2022-06-15
"""
from pathlib import Path
from copy import deepcopy
from string import ascii_uppercase
from collections import defaultdict
from typing import List, Dict, Union

import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb

import enzy_htp.chemical as chem
from enzy_htp.core import _LOGGER, UnsupportedMethod
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
        for mlist in muts.values():
            result.append(mlist.pop())
        return result

    mutations = unique_by_pos(mutations)
    needed: int = n_mutations - len(mutations)
    if needed > 0:
        mutations.extend(
            get_mutations(pdb, needed, mutations, random_state, restrictions)
        )
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

    IMPLEMENTATION: Dict = {"tleap": _mutate_tleap}

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
            raise Exception(
                "Could not generate enough Mutation namedtuple()'s. This is likely due to too many restrictions being present."
            )

        md_keys = list(mut_dict.keys())
        chosen = np.random.choice(md_keys)
        if len(mut_dict[chosen]):
            result.append(np.random.choice(mut_dic[chosen]))
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


def _mutate_tleap(pdb: str, outfile: str, mutations: List[Mutation]) -> None:
    """Underlying implementation of mutation with tleap. Serves as impelementation only, SHOULD NOT
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
    ai = mm.AmberInterface()
    ai.mutate(outfile)
