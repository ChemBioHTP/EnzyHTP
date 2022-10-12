"""Helper submodule that allows for storing of chemical data into a .json file
in enzy_htp/chemical/chemical-database.json. Allows for the python files in this
submodule to be much smaller.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-06-16
"""
import json
import shutil
from pathlib import Path
from typing import Any, Dict

from enzy_htp.core import file_system as fs

_DB_HANDLE: Dict[str, Any] = None
"""Variable for module that holds a dict() with (key, value) pairs of the format 
(variable_name, data), where 'variable_name' is a fully capitalized variable name 
and data is chemical information of any type."""


def load_from_db(var: str) -> Any:
    """Loads the variable value from the 'enzy_htp/chemical/chemical-database.json' database with the
        same key as described by var.

    Args:
        var: The name of the variable whose value should be extracted.

    Returns:
        The associated value for the var or None if the var is not defined.
    """
    global _DB_HANDLE
    if _DB_HANDLE is None:
        db_file: str = str(
            Path(__file__).absolute().parent) + "/chemical-database.json"
        _DB_HANDLE = json.load(open(db_file, "r"))

    return _DB_HANDLE.get(var, None)


def add_to_db(key: str, value: Any) -> None:
    """Saves a new variable value to the 'enzy_htp/chemical/chemical-database.json' databse, ensuring
    that the data is encoded correctly before overwriting the existing database file. NOTE: Does not
        warn if the variable already exists in the databse.

    Args:
        key: The name of the variable to be added.
        var: Data to be saved to the database file.

    Returns:
        Nothing.
    """
    global _DB_HANDLE
    db_file: str = str(
        Path(__file__).absolute().parent) + "/chemical-database.json"
    db_temp: str = str(
        Path(__file__).absolute().parent) + "/chemical-database.json.tmp"
    if _DB_HANDLE is None:
        _DB_HANDLE = json.load(open(db_file, "r"))

    _DB_HANDLE[key] = value

    try:
        json.dump(_DB_HANDLE, open(db_temp, "w"))
        shutil.move(db_temp, db_file)
    except:
        # TODO(CJ): add a warning when it doesnt work.
        pass

    fs.safe_rm(db_temp)
