"""Functions related to decode the target_aa_pattern.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-01-26
"""
import re
from typing import List

from enzy_htp.chemical.residue import (
    ONE_LETTER_CAA_MAPPER,
    ONE_LETTER_AA_MAPPER,
    RESIDUE_VOLUME_MAPPER,
    CAA_CHARGE_MAPPER,
    AA_CHARGE_MAPPER,
    THREE_LETTER_AA_MAPPER,
    RESIDUE_CATEGORIES,
    convert_to_canonical_three_letter,
    convert_to_one_letter,
    convert_to_three_letter,
)
from enzy_htp.core.logger import _LOGGER


def decode_target_aa_pattern(orig_resi: str, pattern: str) -> List[str]:
    """obtain a list of target mutation AA based on a {pattern} and the {orig_resi}
    Args:
        orig_resi: the original residue in the position in its 3-letter code form
        pattern: the pattern that describe the target aa
            syntax:
                1. keywords
                    (current supported keywords)
                    self:       the AA itself
                    all:        all 20 canonical amino acid (AA)
                    larger:     AA that is larger in size according to
                                enzy_htp.chemical.residue.RESIDUE_VOLUME_MAPPER
                    smaller:    AA that is smaller in size
                    similar_size_20: AA that similar is size (cutoff: 20 Ang^3)
                    charge+:    AA that carry more formal positive charge
                    charge-:    AA that carry less formal positive charge
                    charge+1:   AA that carry 1 more positive charge
                    charge-1:   AA that carry 1 less positive charge
                    neutral:    AA that is charge neutral
                    positive:   AA that have positive charge
                    negative:   AA that have negative charge
                    {3-letter}: the AA of the 3-letter name

    Returns:
        a list of 3-letter name of the target AAs"""
    result = []
    keyword_re_pattern = r"\b(?!and|or|not|\(|\))[\d\w+-]+"
    eval_pattern = re.sub(keyword_re_pattern,
                          lambda x: f"set(KEYWORD_DECODER['{x.group(0)}'](orig_resi))",
                          pattern)
    eval_pattern = eval_pattern.replace(" not ", "-")
    eval_pattern = eval_pattern.replace(" or ", "|")
    eval_pattern = eval_pattern.replace(" and ", "&")
    result = eval(eval_pattern)
    return sorted(list(result))


def decode_all(orig_resi: str):
    """decoder for keyword: all"""
    return list(ONE_LETTER_CAA_MAPPER.values())


def decode_larger(orig_resi: str):
    """decoder for keyword: larger"""
    result = filter(
        lambda x: RESIDUE_VOLUME_MAPPER[x] > RESIDUE_VOLUME_MAPPER[convert_to_one_letter(orig_resi)],
        RESIDUE_VOLUME_MAPPER)
    result = map(convert_to_three_letter, result)
    return list(result)


def decode_similar_size_20(orig_resi: str):
    """decoder for keyword: similar_size_20"""
    result = filter(
        lambda x: abs(RESIDUE_VOLUME_MAPPER[x] - RESIDUE_VOLUME_MAPPER[convert_to_one_letter(orig_resi)]) < 20,
        RESIDUE_VOLUME_MAPPER)
    result = map(convert_to_three_letter, result)
    return list(result)

def decode_smaller(orig_resi: str):
    """decoder for keyword: smaller"""
    result = filter(
        lambda x: RESIDUE_VOLUME_MAPPER[x] < RESIDUE_VOLUME_MAPPER[convert_to_one_letter(orig_resi)],
        RESIDUE_VOLUME_MAPPER)
    result = map(convert_to_three_letter, result)
    return list(result)

def decode_charge_p(orig_resi: str):
    """decoder for keyword: charge+"""
    orig_charge = AA_CHARGE_MAPPER[orig_resi]
    result = filter(lambda x: CAA_CHARGE_MAPPER[x] > orig_charge, CAA_CHARGE_MAPPER)
    return list(result)


def decode_charge_n(orig_resi: str):
    """decoder for keyword: charge-"""
    result = filter(lambda x: CAA_CHARGE_MAPPER[x] < AA_CHARGE_MAPPER[orig_resi],
                    CAA_CHARGE_MAPPER)
    return list(result)


def decode_charge_p1(orig_resi: str):
    """decoder for keyword: charge+1"""
    result = filter(lambda x: CAA_CHARGE_MAPPER[x] - AA_CHARGE_MAPPER[orig_resi] == 1,
                    CAA_CHARGE_MAPPER)
    return list(result)


def decode_charge_n1(orig_resi: str):
    """decoder for keyword: charge-1"""
    result = filter(lambda x: CAA_CHARGE_MAPPER[x] - AA_CHARGE_MAPPER[orig_resi] == -1,
                    CAA_CHARGE_MAPPER)
    return list(result)


def decode_neutral(orig_resi: str):
    """decoder for keyword: neutral"""
    return [convert_to_three_letter(i) for i in RESIDUE_CATEGORIES["neutral"]]


def decode_positive(orig_resi: str):
    """decoder for keyword: positive"""
    return [convert_to_three_letter(i) for i in RESIDUE_CATEGORIES["positive"] if i != "H"]


def decode_negative(orig_resi: str):
    """decoder for keyword: negative"""
    return [convert_to_three_letter(i) for i in RESIDUE_CATEGORIES["negative"]]


def decode_self(orig_resi: str):
    """decoder for keyword: self"""
    # canonicalize the 3 letter name (mainly for different protonation state)
    if orig_resi in THREE_LETTER_AA_MAPPER:
        result = convert_to_canonical_three_letter(orig_resi)
    else:
        result = orig_resi

    return [result]


KEYWORD_DECODER = {
    "all": decode_all,
    "larger": decode_larger,
    "smaller": decode_smaller,
    "charge+": decode_charge_p,
    "charge-": decode_charge_n,
    "charge+1": decode_charge_p1,
    "charge-1": decode_charge_n1,
    "neutral": decode_neutral,
    "positive": decode_positive,
    "negative": decode_negative,
    "self": decode_self,
    "similar_size_20": decode_similar_size_20,
    # all direct 3-letter name
    "ALA" : lambda x: ["ALA"],
    "CYS" : lambda x: ["CYS"],
    "ASP" : lambda x: ["ASP"],
    "GLU" : lambda x: ["GLU"],
    "PHE" : lambda x: ["PHE"],
    "GLY" : lambda x: ["GLY"],
    "HIS" : lambda x: ["HIS"],
    "ILE" : lambda x: ["ILE"],
    "LYS" : lambda x: ["LYS"],
    "LEU" : lambda x: ["LEU"],
    "MET" : lambda x: ["MET"],
    "ASN" : lambda x: ["ASN"],
    "PRO" : lambda x: ["PRO"],
    "GLN" : lambda x: ["GLN"],
    "ARG" : lambda x: ["ARG"],
    "SER" : lambda x: ["SER"],
    "THR" : lambda x: ["THR"],
    "SEC" : lambda x: ["SEC"],
    "VAL" : lambda x: ["VAL"],
    "TRP" : lambda x: ["TRP"],
    "TYR" : lambda x: ["TYR"],
    # TODO: add upon need.
}


def check_target_aa_pattern(pattern: str) -> None:
    """check the target aa pattern and give warning as needed
    Current warning:
        - HIS not included"""
    keyword_re_pattern = r"\b(?!and|or|not|\(|\))[\d\w+-]+"
    # give a warning about HIS
    keyword_list = re.findall(keyword_re_pattern, pattern)
    for kw in keyword_list:
        if "charge" in kw:
            _LOGGER.warning(
                f"(Keyword {kw} is used) Note that in keyword decode for target, HIS is not included in the pool due to its subtle titration states."
            )
            return
