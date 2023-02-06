#!/usr/bin/env python

############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import os
import shutil
import sys


def main(args):
    log = logging.getLogger("copy files")
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter("%(message)s"))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    for inputfilename, outputfilename in zip(args[1::2], args[2::2]):
        if os.path.isfile(inputfilename):
            shutil.copyfile(inputfilename, outputfilename)


if __name__ == "__main__":
    main(sys.argv)
