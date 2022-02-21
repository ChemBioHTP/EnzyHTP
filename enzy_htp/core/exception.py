"""TODO: Need a module docstring"""


class MissingEnvironmentElement(Exception):
    """Exception corresponding to a missing executable or environment variable"""
    pass


class InvalidResidueCode(Exception):
	"""Exception corresponding to an invalid one letter or three letter nucleotide being entered"""
	pass
