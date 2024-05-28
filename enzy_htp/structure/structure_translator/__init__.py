"""The structure.structure_translator sub-module revolves around renmaing Atom()'s and Residue()'s 
to a given standard. EnzyHTP assumes the AmberMD naming scheme as a standard but 
the TranslatorBase() represents a template for how to convert to and from AmberMD's naming for an 
arbitrary naming scheme. Users should access this functionality through the translate_structure() 
free function.

Author: Chris Jurich
Date: 2024-02-11
"""

from .translator_base import TranslatorBase

from .api import translate_structure
