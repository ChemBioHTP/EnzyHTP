"""Module holds functions for checking types of arguments. Meant to be used inside of other functions to ensure that
variables are of the correct type. Revolves around check_var_type() function which checks if a supplied variable is 
of a given type.

Author: Chris Jurich, <chris.jurich@vanderbilt.edu>
Date: 2023-06-02
"""
from typing import Any

from .logger import _LOGGER


def _fail_strategy(msg: str, strategy: str) -> None:
    """Implementation method for exit strategy on failure. Takes message and straetgy. Allowed strategies are
    'exit' and 'exception'

    Args:
        msg: The error message as a str().
        strategy: The strategy for exiting. Allowed values are 'exit' and 'exception'.

    Returns:
        Nothing.

    """
    if strategy == 'exit':
        _LOGGER.error(f"{msg} Exiting...")
        exit(1)

    elif strategy == 'exception':
        raise TypeError(err_msg)

    else:
        _LOGGER.error(f"The supplied exit strategy '{strategy}' is not valid. Allowed are 'exit' and 'exception'. Exiting...")
        exit(1)


def check_var_type(var: Any, target_type: Any, fail_strategy: str = 'exit') -> None:
    """Checking if the supplied var is the specified type. Nothing happens if it is correct type. Otherwise will 'exit' or throw
    an exception as specified by the 'fail_strategy' argument.

    Args:
        var: The variable in question.
        target_type: The target type to compare the var to.
        fail_strategy: How type failure is handled. Default is 'exit', the other allowed method is 'exception'.

    Returns:
        Nothing.

    """
    if type(var) == target_type:
        return

    err_msg: str = f"Supplied var '{var}' is NOT a {target_type.__name__}."
    _fail_strategy(err_msg, fail_strategy)
