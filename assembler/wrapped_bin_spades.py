#!/usr/bin/env python

import os
import sys
cur_d = os.path.dirname(os.path.realpath(__file__))

bin = os.path.join(cur_d, "bin/spades")

if len(sys.argv) != 2:
    print "Usage: " + sys.argv[0] + " config file "
os.system (bin  + " " + sys.argv[1] + " " + os.path.join(cur_d, "configs/debruijn/mda_mode.info") + " " +  os.path.join(cur_d, "configs/debruijn/meta_mode.info") + " " +  os.path.join(cur_d, "configs/debruijn/plasmid_mode.info"))
