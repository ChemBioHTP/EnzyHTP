"""Testing the functionality implemented in enzy_htp.chemical.atoms.py

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-17
"""
import enzy_htp.chemical as chem


def is_close(v1: float, v2: float, eps: float = 0.001) -> bool:
    """Checks if two floating point numbers are close enough to be equal.

    Args:
        v1: The first floating point number.
        v2: The second floating point number.
        eps: The allower difference between the two numbers, default is 0.001.

    Returns:
        Whether the two floating point numbers are roughly equivalent.
    """
    return abs(v1 - v2) <= eps


def test_get_h_bond_length_valid_inputs():
    """Testing that the correct h-bond values are supplied for defined values."""
    assert is_close(chem.get_h_bond_length("C"), 1.07)
    assert is_close(chem.get_h_bond_length("N"), 1.0)


def test_get_h_bond_length_invalid_inputs():
    """Testing that the value of -1.0 is supplied for undefined values."""

    for ch in list("ABDEFHGIJKLMOPQRSTUVWXYZ"):
        assert is_close(chem.get_h_bond_length(ch), -1.0)
