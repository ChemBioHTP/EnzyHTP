"""Cluster module includes the support for different clusters.
Add class definations to this module to support more clusters.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Date: 2022-04-13
"""
import os
from importlib import import_module

# import all modules under the dir to support user customization.
__all__ = []
for f in os.listdir(os.path.dirname(__file__)):
    if f != '__init__.py' and f.endswith('.py') and not f.startswith('_'):
        module_name = f.removesuffix('.py')
        __all__.append(module_name)
        import_module('.'+module_name, package='core.clusters')
del f # pylint: disable=undefined-loop-variable

