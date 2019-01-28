#!/usr/bin/env python

import sys

print "#!/bin/sh"
print ""

# vectorized implementations
alg = ["nw", "sg", "sw"]
#stats = ["", "_stats"]
stats = [""]
#par = ["_scan", "_striped", "_diag"]
par = ["_scan"]
#isa = ["_sse2_128_16", "_sse41_128_16", "_avx2_256_16"]
isa = ["_avx2_256_16"]
#blosum = ["40","45","50","62","75","80","90"]
blosum = ["62"]
#threads = ["1", "2", "4", "8", "16", "32", "64", "128", "256"]
threads = ["1", "2", "4", "8", "16", "32", "64"]
for a in alg:
    for s in stats:
        for p in par:
            for i in isa:
                for b in blosum:
                    for t in threads:
                        txt = "OMP_NUM_THREADS="+t+" ./test_openmp -a "+a+s+p+i+" -b blosum"+b+" -f "+sys.argv[1]
                        print "echo '"+txt+"'"
                        print txt

## serial reference implementations
#alg = ["nw", "sg", "sw"]
##stats = ["", "_stats"]
#stats = [""]
##par = ["", "_scan"]
#par = [""]
##blosum = ["40","45","50","62","75","80","90"]
#blosum = ["62"]
#threads = ["1", "2", "4", "8", "16", "32", "64", "128", "256"]
#for a in alg:
#    for s in stats:
#        for p in par:
#            for b in blosum:
#                for t in threads:
#                    txt = "OMP_NUM_THREADS="+t+" ./test_openmp -a "+a+s+p+" -b blosum"+b+" -f "+sys.argv[1]
#                    print "echo '"+txt+"'"
#                    print txt
#
