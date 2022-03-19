"""Defines a logger format using the python logging module. Also creates a singleton _LOGGER for use throughout the module.

Author: Chris Jurich, <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
import os
import logging
import colorlog


def init_logger(dunder_name : str,
                log_outfile: str = None,
                testing_mode: bool = False ,
                start: bool = False) -> logging.Logger:
    """Function for creating the formatted logger. Taken from https://github.com/jyesselm/dreem/dreem/logger.py."""
    log_format = ("[%(asctime)s "
                  "%(name)s "
                  "%(funcName)s] "
                  "%(levelname)s "
                  "%(message)s")
    bold_seq = "\033[1m"
    colorlog_format = f"{bold_seq}" "%(log_color)s" f"{log_format}"
    logger = logging.getLogger(dunder_name)
    # colorlog.basicConfig(format=colorlog_format, datefmt="%H:%M")
    handler = colorlog.StreamHandler()
    handler.setFormatter(
        colorlog.ColoredFormatter(
            colorlog_format,
            datefmt="%H:%M",
            reset=True,
            log_colors={
                "DEBUG": "cyan",
                "WARNING": "yellow",
                "ERROR": "red",
                "CRITICAL": "red,bg_white",
            },
        ))

    logger.addHandler(handler)

    if log_outfile is not None:
        if start:
            if os.path.isfile(log_outfile):
                os.remove(log_outfile)
        fileHandler = logging.FileHandler(log_outfile)
        fileHandler.setFormatter(logging.Formatter(log_format))
        logger.addHandler(fileHandler)

    if testing_mode:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    return logger

  
_LOGGER = init_logger('EnzyHTP', None, start=True)
"""Singleton logging object to log to throught enzy_htp."""
