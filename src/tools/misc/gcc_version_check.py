#!/usr/bin/python

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import os
import time
import subprocess 

gccversion = '4.9.3' if len(sys.argv) < 2 else sys.argv[1] 

cmd_list_to_check = ['gcc', 'g++', 'cc']
excode = 0

for cmd in cmd_list_to_check:
    status, output = subprocess.getstatusoutput(cmd + " -v")
    version = output.split('\n')[-1].split()[2]

    if version != gccversion:
        print("Error: " + cmd + " has version " + version + " instead of " + gccversion)
        excode = 1

exit(excode)

