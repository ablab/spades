
#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys 

def parse_name(name):
      pieces = name.strip().split('_')
      return int(pieces[1])


def find_length(filename):
    inFileName = filename
    inFile = open(inFileName)
    length = 0
    for line in inFile:
        parse = line.strip().split()
        if parse[0] == "10" or parse[0] == "ref0":
            length = max(length, int(parse[3]))
    print "length = " + str(length)
    return length
def find_gaps(filename, length):
    from numpy import ndarray
    a = ndarray((length + 1),int)
    inFileName = filename
    inFile = open(inFileName)
    for line in inFile:
        parse = line.strip().split()
        if parse[0] == "10" or parse[0] == "ref0":
            for j in range(int(parse[1]), int(parse[3])):
                a[j] = 1
    in_gap = 0
    start = 0;
    for i in range(0, length-1):

        if in_gap == 0:
            if a[i] == 0:
                in_gap = 1
                start = i
        if in_gap == 1:
            if a[i] == 1:
                finish = i - 1
                in_gap = 0
                if finish - start > 10:
                    print "gap " + str(start) + " - " + str(finish) +" len: " + str(finish - start)

if len(sys.argv) != 2:
        print "Prints gaps for given .pos file"
        print("Usage: " + sys.argv[0] + " <.pos>")
        sys.exit()
       
length = find_length(sys.argv[1])       
find_gaps(sys.argv[1], length)
#print "*************************************** \nSecond:"
#gaps_after = make_snap(sys.argv[2], length)
#print "####################################### \nCounting gaps ..."
#find__gaps(gaps_before, gaps_after, length)
