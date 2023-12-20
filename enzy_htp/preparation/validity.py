"""Module for the validity check of targets(structure/trajectory/system) of EnzyHTP.
These checks will contain 2 components:
- the target's intrinsic validity (e.g.: whether a structure is chemically reasonable.)
- whether current EnzyHTP can support the chemcial diversity in the target (e.g.: we dont support D/RNA for now.)
The detailed checking logic/criteria may change depending on the version of EnzyHTP.
(i.e.: more diversed systems will be supported)

Science API:
+ is_structure_valid

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2023-12-19
"""
import os
from typing import Tuple

from enzy_htp.core.logger import _LOGGER
from enzy_htp.structure import Structure, StructureEnsemble

def is_structure_valid(stru: Structure, print_report: bool = True) -> Tuple[bool, list]:
    """check 1. the intrinsic validity and 2. EnzyHTP support of
     the chemcial diversity in this structure.
    Args:
        stru: the target Structure()
        print_report: whether print report (reason of failure, how to make right) to logging.
    Returns:
        (boolean of pass/fail, [])
    Note:
        This checker consider the compatiblity of common components of EnzyHTP based workflows
        e.g.: protonate_stru(), mutate_stru(), parameterizer, etc. This mainly serves as a template.
        TODO need to redesign this.
        Concrete workflow author should write their own checker function specialized for the
        workflow."""
    unsupport_reason = [] # (reason, source, suggestion)
    # 1. intrinsic validity
    # TODO develop when found examples.

    # 2. enzyhtp support (TODO add more when encountered)
    future_supp = "this may be supported in future EnzyHTP."
    partial_supp = "this is only partially support. may cause some bug."
    unsupport_mapper = {
        "pdb_io" : {
            "nucleic_acid": partial_supp}, # only partial support, need new class for it.
        "protonate_stru" : {
            "modified_residue": "check https://enzyhtp-doc.readthedocs.io/en/latest/modified_aa.html for manually support this."},
        "mutate_stru": {},
        "AmberParameterizer": {
            "modified_residue" : future_supp, 
            "metalcenters" : "metals are treated as nonbonding from the solvent FF you use for now. We are developing support for MCPB", 
            "nucleic_acid" : future_supp},
        }
    diversity = stru.chemical_diversity
    for source, v in unsupport_mapper.items():
        for reason, suggestion in v.items():
            if reason in diversity:
                unsupport_reason.append((reason, source, suggestion))

    if unsupport_reason:
        if print_report:
            info_str = f"""This structure is not supported for the following reason:
            {'{:<20}'.format('Reason')}{'{:<20}'.format('Source')}{'{:<100}'.format('Suggestion')}{os.linesep}"""
            for reason, source, suggestion in unsupport_reason:
                info_str += f"            {reason:<20}{source:<20}{suggestion:<100}{os.linesep}"
            _LOGGER.warning(info_str)

        return (False, unsupport_reason)
    return (True, [])


def is_stru_ensemble_valid(stru_esm: StructureEnsemble):
    """TODO"""
