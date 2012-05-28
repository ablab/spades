#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import spades_init

spades_init.init()
spades_home = spades_init.spades_home

import quality
quality.main(sys.argv[1:], lib_dir=os.path.join(spades_home, "src/tools/quality/libs"), release_mode=True)
