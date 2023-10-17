"""the I/O submodule for generation, construction, importing, and exporting of various structure related objects:
    + Structure objects: PDBParser
    + Ligand objects: Mol2Parser

Author: Qianzhen (QZ) Shao, <shaoqz@icloud.com>
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-08-01
"""
from .pdb_io import PDBParser
from _interface import StructureParserInterface
from .mol2_io import Mol2Parser
# # import all modules under the dir to support user customization.
# __all__ = []
# for f in os.listdir(os.path.dirname(__file__)):
#     if f != "__init__.py" and f.endswith(".py") and not f.startswith("_"):
#         module_name = f.removesuffix(".py")
#         __all__.append(module_name)
#         import_module("."+module_name, package="core.clusters")
# del f # pylint: disable=undefined-loop-variable
