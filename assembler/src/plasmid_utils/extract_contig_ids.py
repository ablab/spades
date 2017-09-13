#!/usr/bin/python 

import sys
import os
import random
import sets

names = set()
for line in open(sys.argv[1], "r"):
    names.update({int(line.split('_')[1])})
for i in names:
    print i

