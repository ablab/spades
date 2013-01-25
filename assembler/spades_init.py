############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys

spades_home = os.path.abspath(sys.path[0])
spades_version = ''

def init():
    global spades_home
    global spades_version
    global spades_compilation_dir
    if spades_home == "/usr/bin":
        spades_home = "/usr/share/spades"

    sys.path.append(os.path.join(spades_home, "src/spades_pipeline/"))

    spades_version =  open(os.path.join(spades_home, 'VERSION'), 'r').readline()
