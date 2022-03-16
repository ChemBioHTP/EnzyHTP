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
        packages=['enzy_htp',
		'enzy_htp.structure',
		'enzy_htp.chemical',
		'enzy_htp.core',
		],
		install_requirements=requirements
        )
