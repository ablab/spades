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

#filter duplicates workaround(by coverage+length)
if __name__ == "__main__":
#    main_dir = sys.argv[1]
#    dirs = os.listdir(main_dir)
    used = set()
    fullname = join(sys.argv[1])
    to_print = True
    for line in open(fullname, 'r'):
        if line[0] == ">":
            name = line.split()[0]
            lens = name.split('_')[3]
            cov = name.split('_')[5]
            cat = lens+"_"+cov
            if not cat in used and float(cov) > 10 and int(lens) > 10000:
                used.add(cat)
                print line.strip()
                to_print = True
            else:
                to_print = False
        else:
            if to_print:
                print line.strip()
            


