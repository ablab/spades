#!/usr/bin/env python

import sys

print "#!/bin/sh"
print ""

# vectorized implementations
alg = ["nw", "sg", "sw"]
#stats = ["", "_stats"]
stats = [""]
par = ["_scan", "_striped", "_diag"]
#isa = ["_sse2_128_16", "_sse41_128_16", "_avx2_256_16"]
isa = ["_avx2_256_16"]
#blosum = ["40","45","50","62","75","80","90"]
#blosum = ["50", "62"]
#gap = [(10,1), (10,2), (14,2), (40,2)]
blosum_and_gap = [
    ("45",15,2),
    ("50",13,2),
    ("62",11,1),
    ("80",10,1),
    ("90",10,1)]
for b,o,e in blosum_and_gap:
    for a in alg:
        for s in stats:
            for p in par:
                for i in isa:
                    o = str(o)
                    e = str(e)
                    txt = "./tests/test_openmp -a "+a+s+p+i+" -b blosum"+b+" -o "+o+" -e "+e+" -f "+sys.argv[1]+" -i 5"
                    print "echo '"+txt+"'"
                    print txt

