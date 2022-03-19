# TODO documentation
from collections import namedtuple

MutaFlag = namedtuple(
    "MutaFlag", "orig_residue chain_index residue_index target_residue"
)


def mutaflag_to_str(mf: MutaFlag) -> str:
    return f"{mf.orig_residue}{mf.chain_index}{mf.residue_index}{mf.target_residue}"
