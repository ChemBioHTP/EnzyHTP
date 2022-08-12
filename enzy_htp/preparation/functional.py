"""
"""


from pathlib import Path
from typing import Any, List, Union

import numpy as np


from .protonate import protonate_pdb, protonate_missing_elements
from .pdb_line import PDBLine, read_pdb_lines
from enzy_htp.core import file_system as fs

from enzy_htp.structure import (
    compare_structures,
    merge_right,
    structure_from_pdb,
    Ligand,
    ligand_from_pdb,
    Chain,
    Structure,
)


def remove_water(start_pdb: str, end_pdb: str = None) -> str:
    """ """
    pdb_lines: List[PDBLine] = read_pdb_lines(start_pdb)
    pdb_lines = list(
        filter(lambda pl: not (pl.is_water() or pl.is_CRYST1()), pdb_lines)
    )
    mask: List[bool] = [True] * len(pdb_lines)

    for pidx, pl in enumerate(pdb_lines[:-1]):
        if pl.is_TER() and pdb_lines[pidx + 1].is_TER():
            mask[pidx + 1] = False

    pdb_lines = np.array(pdb_lines)[mask]

    outfile: str = start_pdb
    if end_pdb:
        outfile = end_pdb

    fs.write_lines(outfile, list(map(str, pdb_lines)))
    # TODO(CJ): may need to reset line numbers here before writing.
    return outfile


def protonate(
    start_pdb: str, end_pdb: str = None, ph: float = 7.0, ffout: str = "AMBER"
) -> str:
    # steps
    # 1. protonate the .pdb file
    # 2. protonate non  elements
    # 3. merge them back in
    files_to_delete: List[str] = list()
    pqr = str(Path(start_pdb).with_suffix(".pqr.pdb"))
    pqr_log = str(Path(start_pdb).with_suffix(".pqr.log"))
    files_to_delete.append(pqr)
    files_to_delete.append(pqr_log)
    protonate_pdb(start_pdb, pqr)
    new_structure: Structure = protonate_missing_elements(
        start_pdb, pqr, str(Path(start_pdb).parent)
    )

    outfile: str = start_pdb
    if end_pdb:
        outfile = end_pdb

    new_structure.to_pdb(outfile)

    _ = list(map(fs.safe_rm, files_to_delete))

    return outfile


def prepare_from_pdb(start_pdb: str, outdir: str, ph: float = 7.0) -> str:
    """
    """
    fs.safe_mkdir(outdir)
    outfile: str = f"{outdir}/{Path(start_pdb).stem}_prepared.pdb"
    remove_water(start_pdb, outfile)
    protonate(outfile, outfile, ph)

    return outfile
