"""the I/O submodule for generation/construction of Structure objects from files and other data structures (?should it be here) and exporting it vice versa

Author: shaoqz, <shaoqz@icloud.com>
Date: 2022-08-01
"""
from .pdb_io import PDBParser
# # import all modules under the dir to support user customization.
# __all__ = []
# for f in os.listdir(os.path.dirname(__file__)):
#     if f != '__init__.py' and f.endswith('.py') and not f.startswith('_'):
#         module_name = f.removesuffix('.py')
#         __all__.append(module_name)
#         import_module('.'+module_name, package='core.clusters')
# del f # pylint: disable=undefined-loop-variable