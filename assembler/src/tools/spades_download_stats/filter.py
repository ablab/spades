#!/usr/bin/python

#***************************************************************************
##* Copyright (c) 2011-2014 Saint-Petersburg Academic University
##* All Rights Reserved
##* See file LICENSE for details.
#****************************************************************************


import sys

f = open(sys.argv[1])

l = []
for line in f:
    l.append(line.split()[0])
    

s = set(l)
print(str(len(l)) + " (" + str(len(s)) +  ")")



