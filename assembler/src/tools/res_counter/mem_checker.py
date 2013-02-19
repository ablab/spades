#!/usr/bin/python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import os
import time
import commands


maxmem = 0

while True:
    status, output = commands.getstatusoutput("du -h --max-depth=0")
    memstr = output.split()[0]

    if memstr[-1] == "K":
        mem = float(memstr[0:-1:1]) * 1024
    elif memstr[-1] == "M":
        mem = float(memstr[0:-1:1]) * 1024 * 1024
    elif memstr[-1] == "G":
        mem = float(memstr[0:-1:1]) * 1024 * 1024 * 1024
    elif memstr[-1] == "T":
        mem = float(memstr[0:-1:1]) * 1024 * 1024 * 1024 * 1024
    else:
        mem = float(memstr)

    if mem > maxmem:
        maxmem = mem

    outf = open("maxmem.txt", "w")
    outf.write(str(maxmem))
    outf.close()

    time.sleep(10)
    
