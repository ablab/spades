#!/usr/bin/env python
from __future__ import print_function
try:
    from itertools import izip as zip
except ImportError:
    pass
import os
import sys

print(" ".join(sys.argv))
if len(sys.argv) < 3:
    print("Usage: %s <output file> (<contig .id file>)+" % sys.argv[0])
    print(".mpl file should be at the same folders with same names")
    sys.exit(1)


with open(sys.argv[1], "w") as out:
    files = sys.argv[2:]

    for f in files:
        if (f[-3:] != ".id"):
            print(".id files are supposed as input")
            sys.exit(2)

        print("Processing abundances from %s" % f)
        name = os.path.splitext(os.path.basename(f))[0]

        with open(f, "r") as ctg_id, open(f[:-3] + ".mpl", "r") as ctg_mpl:
            for cid, cmpl in zip(ctg_id, ctg_mpl):
                print(name + "-" + cid.strip() + " " + cmpl.strip(), file=out)
