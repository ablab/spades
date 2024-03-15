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
    contigs = sorted(glob.glob(os.path.join(contig_dir, "*.fast?")) + glob.glob(os.path.join(contig_dir, "*.gfa")) + glob.glob(os.path.join(contig_dir, "*.paths")))
    except_list = []
    for c in contigs:
        if c.find("latest_") != -1:
            print("Will not remove symlink " + c)
            except_list.append(c)
            print("Will not remove " + os.path.realpath(c))
            except_list.append(os.path.realpath(c))

    for c in contigs:
        if c in except_list:
            print("Skipping " + c)
            continue
        diff = datetime.datetime.now() - datetime.datetime.fromtimestamp(os.path.getmtime(c))
        if diff.days >= days_difference:
            try:
                os.remove(c)
                print("Removed " + c)
            except OSError:
                print("Error occured while trying to remove " + c)



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

            
