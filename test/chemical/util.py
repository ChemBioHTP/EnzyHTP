"""Testing utilities for the testing of enzy_htp.chemical 

Author: Chris Jurich <chris.jurich@vanderbilt.edu
Date: 2022-03-19
"""

from typing import List

def all_caps( raw_vals : List[str] ) -> bool:
    """Helper method that checks all str() in a list() are capitalized."""
    for rv in raw_vals:
        if not rv.upper():
            return False 
    return True

