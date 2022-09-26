"""Module holds functions for math calculations. 

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2022-09-26
"""
from .logger import _LOGGER

def check_valid_ph(ph: float) -> None:
    """Helper function that checks the pH is on the range: [0.00, 14.00]. Give a warning if not."""
    if ph < 0 or ph > 14:
        _LOGGER.warning(f"assigned pH: {ph:.2f} out of range: [0.00,14.00]")
