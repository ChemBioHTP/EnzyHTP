"""
core module for enzy_htp. Responsibilities include file system management, environment management, exceptions and logging.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
from .logger import _LOGGER
from .env_manager import EnvironmentManager
from .exception import (
    MissingEnvironmentElement,
    InvalidResidueCode,
    UnsupportedFileType,
    UnsupportedMethod,
    InvalidPH,
    InvalidMutationRestriction,
)
from .file_system import (
    safe_rm,
    safe_rmdir,
    safe_mkdir,
    base_file_name,
    get_current_time,
    lines_from_file,
    write_lines,
    get_file_ext,
    write_data,
)
