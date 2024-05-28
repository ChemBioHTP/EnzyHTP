"""Free function for translating the names of a Structure()'s Residue() and Atom() objects. Should be called
by users, supports the conversion of names between EnzyHTP standard (AmberMD) and:

    + Rosetta

as well as converting back to EnzyHTP standard. All functionality is handled through `translate_structure()`.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2024-02-11
"""
from typing import Dict, Union

from enzy_htp.core import _LOGGER

from ..atom import Atom
from ..residue import Residue
from ..structure import Structure

from .translator_base import TranslatorBase
from .rosetta_translator import RosettaTranslator


def translate_structure(stru:Structure, start_naming:str=None, end_naming:str=None) -> None:
    """Translates names of the supplied Structure()'s Residue() and Atom() objects to EnzyHTP standard or 
    a supported naming convention. Exactly one of start_naming or end_naming scheme must be supplied. It is assumed that
    the other, blank parameter is EnzyHTP standard/AmberMD. end_naming should be left blank if a Structure() is to be 
    converted to EnzyHTP standard and start_naming should be left blank if a Structure() is already in EnzyHTP standard naming.
    Function performs basic checks on supplied naming parameters and will error if they are invalid. Naming is always done
    in place.

    Args:
        stru: The Structure() to rename.
        start_naming: The starting naming scheme as a str().
        end_naming:  The ending naming scheme as a str().

    Returns:
        Nothing.
    """
    
    err_msg = str()
    if start_naming is None and end_naming is None:
        err_msg = f"To translate a Structure() you must supply either start_naming or end_naming! Both are empty!"

    elif start_naming is not None and end_naming is not None:
        err_msg = f"To translate a Structure() you must supply either start_naming or end_naming! Both have non-empty values!"
    
    else:
        if start_naming:
            translator = TRANSLATORS.get(start_naming, None)
            if translator is None:
                err_msg = f"The supplied start_naming scheme {start_naming} is not supported. Allowed values are {', '.join(TRANSLATORS.keys())}"
            else:
                translator.to_standard(stru)
    
        if end_naming:
            translator = TRANSLATORS.get(end_naming, None)
            if translator is None:
                err_msg = f"The supplied end_naming scheme {end_naming} is not supported. Allowed values are {', '.join(TRANSLATORS.keys())}"
            else:
                translator.from_standard(stru)
    
    if err_msg:
        _LOGGER.error(err_msg)
        raise ValueError(err_msg) 


TRANSLATORS:Dict[str, TranslatorBase] = {
    'rosetta': RosettaTranslator()
}
"""Mapper that holds instances of all the supported Translators."""
