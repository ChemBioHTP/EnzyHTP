# python imports
import sys 
import logging


# enzy_htp imports
from .structure import Structure
from .core import check_compatability

def init_logging() -> None:
	logging.basicConfig(level=logging.INFO, handlers=[logging.StreamHandler(sys.stdout)])
	

init_logging()
check_compatability()
