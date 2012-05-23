############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys

spades_home = os.path.abspath(sys.path[0])
spades_version = ''
spades_build_dir = os.path.join(spades_home, "build")

def init():
    global spades_home
    global spades_version
    global spades_compilation_dir
    if spades_home == "/usr/bin":
        spades_home = "/usr/share/spades"
        spades_build_dir = os.path.join(os.getenv('HOME'), '.spades')

    sys.path.append(os.path.join(spades_home, "src/tools/spades_pipeline/"))
    sys.path.append(os.path.join(spades_home, "src/tools/quality/"))
    sys.path.append(os.path.join(spades_home, "src/tools/quality/libs"))

    spades_version =  open(os.path.join(spades_home, 'VERSION'), 'r').readline()
