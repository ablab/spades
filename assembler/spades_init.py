#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
from os.path import abspath, dirname, realpath, join, isfile

source_dirs = ["", "truspades", "common"]

# developers configuration
spades_home = abspath(dirname(realpath(__file__)))
bin_home = join(spades_home, 'bin')
python_modules_home = join(spades_home, 'src')
ext_python_modules_home = join(spades_home, 'ext', 'src', 'python_libs')
spades_version = ''


def init():
    global spades_home
    global bin_home
    global python_modules_home
    global spades_version
    global ext_python_modules_home

    # users configuration (spades_init.py and spades binary are in the same directory)
    if isfile(os.path.join(spades_home, 'spades')):
        install_prefix = dirname(spades_home)
        bin_home = join(install_prefix, 'bin')
        spades_home = join(install_prefix, 'share', 'spades')
        python_modules_home = spades_home
        ext_python_modules_home = spades_home

    for dir in source_dirs:
        sys.path.append(join(python_modules_home, 'spades_pipeline', dir))

    spades_version = open(join(spades_home, 'VERSION'), 'r').readline().strip()


if __name__ == '__main__':
    spades_py_path = join(dirname(realpath(__file__)), 'spades.py')
    sys.stderr.write('Please use ' + spades_py_path + ' for running SPAdes genome assembler\n')