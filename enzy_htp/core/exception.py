"""Defines enzy_htp exceptions that should be used instead of generic python exceptions.

Author: Chris Jurich, <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""


class MissingEnvironmentElement(Exception):
    """Exception corresponding to a missing executable or environment variable."""

    pass


class InvalidResidueCode(Exception):
    """Exception corresponding to an invalid one letter or three letter nucleotide being entered."""

    pass


class UnsupportedFileType(Exception):
    """Exception corresponding to a file type that enzy_htp does not support."""

    pass


class UnsupportedMethod(Exception):
    """Exception corresponding to a method that is NOT implemented and/or supported."""

    pass


class InvalidMutationRestriction(Exception):
    """Exception corresponding to an invalid mutation restriction in enzy_htp.mutation.mutation_restrictions.py."""

    pass


class ResidueDontHaveAtom(Exception):
    """Exception corresponding to asking a residue for an atom that it doesn't contain
    takes the residue and the query atom name to constuct"""

    def __init__(self, residue, atom_name, *args) -> None:
        super().__init__(*args)
        self.residue = residue
        self.atom_name = atom_name


class InvalidMutationPatternSyntax(Exception):
    """Exception corresponding to an invalid mutation syntax in enzy_htp.mutation.mutation_pattern"""
    pass


class InvalidMutationFlagSyntax(InvalidMutationPatternSyntax):
    """Exception corresponding to an invalid mutation syntax in enzy_htp.mutation.mutation"""
    pass


class InvalidMutation(Exception):
    """Exception corresponding to an invalid mutation in enzy_htp.mutation.mutation"""
    pass


class tLEaPError(Exception):
    """Exception corresponding runtime error of tleap in enzy_htp._interface.amber_interface.run_tleap()
    contains a list of error information"""

    def __init__(self, error_info_list, *args) -> None:
        super().__init__(*args)
        self.error_info_list = error_info_list

    @property
    def error_info_str(self) -> str:
        return "\n".join(self.error_info_list)


class IndexMappingError(Exception):
    """Exception corresponding index mapping error for residue in enzymes or atoms in molecules"""
    pass


class BadMutantStructure(Exception):
    """Exception corresponding to structure that didn't pass the check in
    enzy_htp.mutation.api.check_mutant_stru"""
    pass

class InconsistentMDEngine(Exception):
    """Exception corresponding to inconsistency of engine between MD building
    blocks in enzy_htp.geometry.sampling"""
    pass

class FileFormatError(Exception):
    """Exception corresponding to wrong file format. Mainly used in file parsing in
    StructureParserInterface concrete classes"""
    pass

class AmberMDError(Exception):
    """Exception corresponding runtime error of sander/pmemd in enzy_htp._interface.amber_interface.AmberMDStep
    contains a list of error information"""

    def __init__(self, error_info_list, *args) -> None:
        super().__init__(*args)
        self.error_info_list = error_info_list

    @property
    def error_info_str(self) -> str:
        return "\n".join(self.error_info_list)

class WrongTopology(Exception):
    """Exception corresponding to wrong topology. Mainly used in StructureConstraint"""
    pass

class GaussianError(Exception):
    """Exception corresponding runtime error of g16 in enzy_htp._interface.gaussian_interface
    contains a list of error information"""

    def __init__(self, error_info_list, *args) -> None:
        super().__init__(*args)
        self.error_info_list = error_info_list

    @property
    def error_info_str(self) -> str:
        return "\n".join(self.error_info_list)

class RosettaError(Exception):
    """Exception corresponding runtime error of Rosetta in enzy_htp._interface.rosetta_interface
    contains a list of error information"""

    def __init__(self, error_info_list, *args) -> None:
        super().__init__(*args)
        self.error_info_list = error_info_list

    @property
    def error_info_str(self) -> str:
        return "\n".join(self.error_info_list)

class ShrapnelChildError(Exception):
    """Exception corresponding to incomplete shrapnel children runs"""
    pass

class AddPDBError(Exception):
    """Exception corresponding to incomplete prmtop file due to no add_pdb"""
    pass
