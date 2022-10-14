#!/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

with open('requirements.txt', 'r') as f:
	requirements = f.read().splitlines()

#TODO add entry points for executables
setup(name='enzy_htp',
        version='0.1.1',
        description='TODO',
        author='EnzyHTP Authors',
        author_email='zhongyue.yang@vanderbilt.edu',
        include_package_data=True,
        packages=['enzy_htp',
		'enzy_htp._config',
		'enzy_htp._interface',
		'enzy_htp.core',
		'enzy_htp.chemical',
		'enzy_htp.structure',
		'enzy_htp.structure.structure_io',
		'enzy_htp.structure.structure_operation',
		'enzy_htp.preparation',
		'enzy_htp.mutation',
		'enzy_htp.molecular_mechanics',
		'enzy_htp.quantum_mechanics',
		],
		install_requires=requirements,
        )
