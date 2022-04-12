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
        version='0.1.0',
        description='TODO',
        author='EnzyHTP Authors',
        author_email='zhongyue.yang@vanderbilt.edu',
        packages=['enzy_htp'],
		py_modules=['enzy_htp.core',
		'enzy_htp.chemical',
		'enzy_htp.structure',
		'enzy_htp.preparation',
		'enzy_htp.molecular_mechanics',
		'enzy_htp.quantum_mechanics',
		'enzy_htp.data',
		],
		install_requirements=requirements,
        include_package_data=True
        )
