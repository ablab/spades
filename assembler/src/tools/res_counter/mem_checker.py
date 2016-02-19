#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import os
import time
import commands
from os.path import isfile

stopper_fname = 'mem_checker_stopper.txt'
maxmem = 0
timeout = 5
dir_to_check = '.' if len(sys.argv) < 2 else sys.argv[1] 

if not isfile(stopper_fname):
  open(stopper_fname, 'w').close()

while isfile(stopper_fname):
    status, output = commands.getstatusoutput("du -s " + dir_to_check)
    memstr = output.split('\n')[-1].split()[0]

    mem = 0
    try:    
        mem = float(memstr) * 1024
    except ValueError:
	pass

    if mem > maxmem:
        maxmem = mem

    with open("maxmem.txt", "w") as f:
        f.write("%db %dM %.1fG\n" % (maxmem // 1,  maxmem // 1024 // 1024, maxmem // 1024 // 1024 / 1024.0))

    time.sleep(timeout)
