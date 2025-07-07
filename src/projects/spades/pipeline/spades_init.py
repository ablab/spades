#!/usr/bin/env python3

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
from os.path import abspath, dirname, realpath, join, isfile

source_dirs = ["", "stages", "common", "executors", "scripts"]

# developers configuration
spades_home = abspath(dirname(realpath(__file__)))
spades_root = abspath(join(spades_home, "../../../../"))
config_dirs = [(join(spades_root, "src", "projects", "spades", "configs"), "debruijn"),
               (join(spades_root, "src", "projects", "ionhammer", "configs"), "ionhammer"),
               (join(spades_root, "src", "projects", "hammer", "configs"), "hammer"),
               (join(spades_root, "src", "projects", "corrector", "configs"), "corrector")]
bin_home = join(spades_root, "bin")
python_modules_home = join(spades_root, "src", "projects", "spades", "pipeline")
ext_python_modules_home = join(spades_root, "ext", "src", "python_libs")
spades_version = ""


def init():
    global spades_home
    global bin_home
    global config_dirs
    global python_modules_home
    global spades_version
    global ext_python_modules_home
    global spades_root

    if isfile(os.path.join(spades_home, "spades-core")):
        # users configuration (spades_init.py and spades binary are in the same directory)
        install_prefix = dirname(spades_home)
        bin_home = join(install_prefix, "bin")
        spades_home = join(install_prefix, "share", "spades")
        config_dirs = [(join(spades_home, "configs"), "./")]
        python_modules_home = spades_home
        ext_python_modules_home = spades_home
        spades_version = open(join(spades_home, "VERSION"), 'r').readline().strip()
    else:
        spades_version = open(join(spades_root, "VERSION"), 'r').readline().strip()

    sys.path.append(join(python_modules_home, "spades_pipeline"))
    for dir in source_dirs:
        sys.path.append(join(python_modules_home, "spades_pipeline", dir))


if __name__ == "__main__":
    spades_py_path = join(dirname(realpath(__file__)), "spades.py")
    sys.stderr.write("Please use " + spades_py_path + " for running SPAdes genome assembler\n")
