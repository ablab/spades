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

gccversion = '4.9.3' if len(sys.argv) < 2 else sys.argv[1] 

cmd_list_to_check = ['gcc', 'g++', 'cc']
excode = 0

for cmd in cmd_list_to_check:
    status, output = commands.getstatusoutput(cmd + " -v")
    version = output.split('\n')[-1].split()[2]

    if version != gccversion:
        print("Error: " + cmd + " has version " + version + " instead of " + gccversion)
        excode = 1

exit(excode)

