#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

#Remove contigs older than given number of days
# options: <dir> <days>

import sys
import os
import time
import glob
import datetime


def remove_in_dir(contig_dir, days_difference):
    contigs = sorted(glob.glob(os.path.join(contig_dir, "*.fasta")))
    for c in contigs:
        diff = datetime.datetime.now() - datetime.datetime.fromtimestamp(os.path.getmtime(c))
        if diff.days >= days_difference:
            os.remove(c)
            print("Removed " + c)


def remove_recursive(storage_dir, days_difference):
    remove_in_dir(storage_dir, days_difference)
    files = sorted(glob.glob(os.path.join(storage_dir, "*")))
    for f in files:
        if os.path.isdir(f):
            remove_recursive(f, days_difference)


if len(sys.argv) != 3:
    print("Run: " + sys.argv[0] + " <dir> <days> ")
    exit(1)

remove_recursive(sys.argv[1], int(sys.argv[2]))

            
