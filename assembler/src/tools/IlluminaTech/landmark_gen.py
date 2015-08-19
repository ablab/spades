############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import random
import sys

__author__ = 'anton'

def MutateReference(s, p):
    list = [c for c in s]
    for i in range(len(list)):
        if list[i] == 'C' and random.random() < p:
            list[i] = random.choice("AGT")
    return "".join(list)

def run(input, output, p):
    lines = open(input, "r").readlines()
    mutated = MutateReference("".join([line.strip() for line in lines[1:]]), p)
    handle = open(output, "w");
    handle.write(lines[0])
    handle.write(mutated)
    handle.close()

run(sys.argv[1], sys.argv[2], float(sys.argv[3]))