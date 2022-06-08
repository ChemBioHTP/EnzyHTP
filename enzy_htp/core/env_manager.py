"""Defines an EnvironmentManager class that checks for applications and environment variables that enzy_htp needs. Also stores values for these variables 
and serves as interface for running external applications.

Author: Chris Jurich, <chris.jurich@vanderbilt.edu>
Date: 2022-02-12
"""
import os
import shutil
from typing import List
from subprocess import run

from .logger import _LOGGER
from .exception import MissingEnvironmentElement


class EnvironmentManager:
    """Serves as general interface between module and the current computer environment.
    Checks whether given applications and environment variables are set in the current environment.
    After check, stores names of executables.
    Serves as interfrace for running commands on system.


    Attributes:
            env_vars_: A list of strings containing environment variables to check for.
            executables_: a list of strings containing executables to check for.
            missing_env_vars_: A list of strings corresponding to environment variables that are missing.
            missing_executables_: A list of strings corresponding to executables that are missing.
    """

    def __init__(self, **kwargs):
        """Initializes object, optionally with starting environment variables and executables."""
        self.env_vars_ = kwargs.get("env_vars", [])
        self.executables_ = kwargs.get("executables", [])
        self.mapper = dict()
        self.missing_env_vars_ = []
        self.missing_executables_ = []

    def add_executable(self, exe_name: str) -> None:
        """Adds the name of an executable to check for."""
        self.executables_.append(exe_name)

    def add_env_var(self, env_var: str) -> None:
        """Adds the name of an environment variable to check for."""
        self.env_vars_.append(env_var)

    def check_env_vars(self) -> None:
        """Checks which environment variables are defined, storing those which are not defined."""
        for env_var in self.env_vars_:
            if os.getenv(env_var) is None:
                self.missing_env_vars_.append(env_var)

    def __exe_exists(self, exe_name: str) -> bool:
        """Helper method that checks if executable exists in current environment."""
        full_path = os.path.expandvars(exe_name)
        return shutil.which(exe_name) is not None

    def check_executables(self) -> None:
        """Checks which executables are available in the system, storing paths to those which exist or noting if they are not found."""
        for exe in self.executables_:
            if not self.__exe_exists(exe):
                self.missing_executables_.append(exe)
            else:
                self.mapper[exe] = shutil.which(exe)

    def display_missing(self) -> None:
        """Displays a list of missing environment variables and exectuables to the logger. Should be called after .check_environment() and .check_env_vars()."""
        if not self.is_missing():
            _LOGGER.info("Environment has all required elements!")
            return
        _LOGGER.warning("Environment is missing some required elements...")

        if len(self.missing_executables_):
            _LOGGER.warning("\tMissing excecutables:")
            for me in self.missing_executables_:
                _LOGGER.warning(f"\t\t{me}")

        if len(self.missing_env_vars_):
            _LOGGER.warning("\tMissing environment variables:")
            for mev in self.missing_env_vars_:
                _LOGGER.warning(f"\t\t{mev}")

    def check_environment(self) -> None:
        """Preferred client method for validating environment. Performs checks and logs output."""
        _LOGGER.info("Checking environment for required elements...")
        self.check_env_vars()
        self.check_executables()
        self.display_missing()
        _LOGGER.info("Environment check completed!")

    def reset(self) -> None:
        """Resets internal lists of env vars and executables."""
        self.executables_ = []
        self.missing_executables_ = []
        self.env_vars_ = []
        self.missing_env_vars_ = []

    def is_missing(self) -> bool:
        """Checks if any executables or environment variables are missing."""
        return len(self.missing_executables_) or len(self.missing_env_vars_)

    def run_command(self, exe: str, args: List[str]) -> List[str]:
        """Interface to run a command with the exectuables specified by exe as well as a list of arguments."""
        cmd = f"{self.mapper.get(exe,exe)} {' '.join(args)}"
        if exe in self.missing_executables_ or not self.__exe_exists(exe):
            _LOGGER.error(
                f"This environment is missing '{exe}' and cannot run the command '{cmd}'"
            )
            _LOGGER.error(f"Exiting...")
            exit(1)
        _LOGGER.info(f"Running command: '{cmd}'...")
        try:
            result = run(cmd, shell=True, capture_output=True)
            res_lines = list(
                map(lambda ss: ss.decode("utf-8"), result.stdout.splitlines())
            )
            _LOGGER.info(f"Command run!")
        except Exception as e:
            _LOGGER.error(
                f"Following error was raised during command excecution: '{e}'. Exiting..."
            )
            exit(1)
        return res_lines

    def __getattr__(self, key: str) -> str:
        """Allows accession into acquired executables."""
        if key not in self.mapper and key in self.executables_:
            _LOGGER.error(
                f"Executable '{key}' is in list of executables to check but has not been searched for yet. Call .check_environment() first. Exiting..."
            )
            exit(1)
        return self.mapper[key]
