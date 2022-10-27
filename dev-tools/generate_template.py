#!/usr/bin/env python3
import os
import argparse
from datetime import datetime


def setup_params():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file_name',
                        type=str,
                        required=True,
                        help='name of python file to be created')
    parser.add_argument(
        '--is_class',
        action='store_true',
        help='Flag to indicate if file contains a class. False by default.')
    parser.add_argument('--description',
                        required=True,
                        type=str,
                        help='description string for the module.')
    return parser.parse_args()


def get_dev_name() -> str:
    result = os.popen('git config user.name').read().splitlines()
    if not result:
        raise Exception(
            "Developer's name is not set. Set with 'git config user.name \"first_name last_name\""
        )
    return result[0]


def get_dev_email():
    result = os.popen('git config user.email').read().splitlines()
    if not result:
        raise Exception(
            "Developer's email is not set. Set with 'git config user.email \"<uname>@email.com\""
        )
    return result[0]


def prepare_name(rawfname: str) -> str:
    if rawfname.endswith('.py'):
        return rawfname
    else:
        return rawfname + '.py'


def class_from_fname(fname):
    fname = fname.replace('.py', '')
    tks = fname.split('_')
    result = str()
    for tk in tks:
        tk = list(tk)
        tk[0] = tk[0].upper()
        result += ''.join(tk)

    return result


def main(params):

    ESCAPED_TRIPLE = "\"\"\""
    TODAY = datetime.now().strftime("%Y-%m-%d")
    fname = prepare_name(params.file_name)
    fh = open(fname, 'w')
    fh.write(f"""{ESCAPED_TRIPLE}{params.description}

Author: {get_dev_name()}, <{get_dev_email()}>
Date: {TODAY}
{ESCAPED_TRIPLE}

""")

    if params.is_class:
        cname = class_from_fname(fname)

        fh.write(f"class {cname}:\n")
        fh.write(f"\tpass\n")

    fh.close()


if __name__ == '__main__':
    main(setup_params())
