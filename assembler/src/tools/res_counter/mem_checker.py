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


maxmem = 0

while True:
    status, output = commands.getstatusoutput("du --max-depth=0")
    memstr = output.split('\n')[-1].split()[0]

    mem = 0
    try:    
        mem = float(memstr) * 1024
    except ValueError:
	pass

    if mem > maxmem:
        maxmem = mem

    outf = open("maxmem.txt", "w")
    outf.write(str(maxmem))
    outf.close()

    time.sleep(10)
    
