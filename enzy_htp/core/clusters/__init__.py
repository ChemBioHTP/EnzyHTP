"""Cluster module includes the support for different clusters.
Add class definations to this module to support more clusters.

Author: Qianzhen (QZ) Shao <qianzhen.shao@vanderbilt.edu>
Date: 2022-04-13
"""
import os
from importlib import import_module

from ._interface import ClusterInterface

# import all modules under the dir to support user customization.
__all__ = []
for f in os.listdir(os.path.dirname(__file__)):
    if not f.startswith('_') and f.endswith('.py'):
        module_name = f.removesuffix('.py')
        __all__.append(module_name)
        import_module('.'+module_name, package='enzy_htp.core.clusters')
del f # pylint: disable=undefined-loop-variable
