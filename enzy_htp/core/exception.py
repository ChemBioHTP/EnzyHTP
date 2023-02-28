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
