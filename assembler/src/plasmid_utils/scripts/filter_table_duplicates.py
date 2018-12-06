#!/usr/bin/python
import sys
import os
import string
import re
import subprocess
import datetime
import fastaparser
from os.path import join
from genericpath import isdir, exists

min_l = 1000

if __name__ == "__main__":
#    main_dir = sys.argv[1]
#    dirs = os.listdir(main_dir)
    used = set()
    fullname = join(sys.argv[1])
    for line in open(fullname, 'r'):
        if len(line) < 10:
            continue
        name = line.split()[0]
        lens = name.split('_')[3]
        cov = name.split('_')[5]
        cat = lens+"_"+cov
        if int(lens) < 500:
            continue
        if not cat in used:
            used.add(cat)
            print line.strip()
            


