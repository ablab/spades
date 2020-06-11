#!/usr/bin/env python

############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import os
from site import addsitedir

import spades_init
spades_init.init()
spades_home = spades_init.spades_home
ext_python_modules_home = spades_init.ext_python_modules_home

addsitedir(ext_python_modules_home)

import options_storage
from spades import main as spades_main

if __name__ == '__main__':
    args = sys.argv
    if len(args) > 1:
        args += ["--rnaviral", "--custom-hmms", os.path.join(spades_home, options_storage.coronaspades_hmms)]
    spades_main(args)
