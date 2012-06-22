#!/usr/bin/env python

import os
import shutil
import sys
import getopt
import glob

import spades_init

spades_init.init()
spades_home = spades_init.spades_home

import support
from process_cfg import *

def prepare_config(filename, K):
    subst_dict = dict()
    subst_dict["K"] = K
    substitute_params(filename, subst_dict)


if len(sys.argv) != 2:
    print ("Usage: " + sys.argv[0] + " K")
    exit(0)

config = "configs/debruijn/config.info"
K = sys.argv[1]
print ("Changing to K = " + K)
prepare_config(config, K)
    


