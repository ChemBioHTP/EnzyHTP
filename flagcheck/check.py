import os
import shutil
import json
import pathlib
import pickle
import random
import sys
import time
from absl import app
from absl import flags
from absl import logging

flags.DEFINE_string('amber_binary_path', shutil.which('am1bcc'),   'Path to the Amber executable.')
flags.DEFINE_string('gaussian_binary_path', shutil.which('g09'),  'Path to the Gaussian executable.')
flags.DEFINE_string('tleap_binary_path', shutil.which('tleap'),   'Path to the tleap executable.')

FLAGS = flags.FLAGS


#def _check_flag(flag_name: str, other_flag_name: 
#                str,should_be_set: bool):
#    if should_be_set != bool(FLAGS[flag_name].value):
#        verb = 'be' if should_be_set else 'not be'
#        raise ValueError(f'{flag_name} must {verb} set when running with '                
#                         f'"--{other_flag_name}={FLAGS[other_flag_name].value}".')


def main(argv):
        if len(argv) > 1:
            raise app.UsageError('Too many command-line arguments')

        for tool_name in ('amber','tleap','gaussian'):
            if not FLAGS[f'{tool_name}_binary_path'].value:
                raise ValueError(f'Could not find path to the "{tool_name}"'
                                 'binary. Make sure it is installed on your system')
            else: 
                print("software installed")

#def is_tools_available(name):
#    return shutil.which(name) is not None


"""
##TODO: requirement
if __name__ == '__main__':
    flag.mark_flags_as_required([
        'g09', 
        'am1bcc',
    ])
"""


app.run(main)
