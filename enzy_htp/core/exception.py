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


class InvalidPH(Exception):
    """Exception corresponding to an invalid pH not in the range 0-14."""

    pass
