#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

#Force symlinks to master contigs in storage
# options: <dir>
#RUN ONLY WHEN ALL TESTS ARE GREEN

import sys
import os
import glob

def force_symlinks_for_suffix(contig_dir, suffix, name):
    contigs = sorted(glob.glob(os.path.join(contig_dir, "*" + suffix + "_master.fasta")))
    if len(contigs) == 0:
        return
    clink = os.path.join(contig_dir, name)
    if os.path.islink(clink):
        os.remove(clink)
    os.symlink(os.path.basename(contigs[-1]), clink)
    print("Forced symlink to " + contigs[-1] + " -> " + clink)


def force_symlinks(contig_dir):
    force_symlinks_for_suffix(contig_dir, "ctg", "latest_contigs.fasta")
    force_symlinks_for_suffix(contig_dir, "scafs", "latest_scaffolds.fasta")


def force_symlinks_recursive(storage_dir):
    force_symlinks(storage_dir)
    files = sorted(glob.glob(os.path.join(storage_dir, "*")))
    for f in files:
        if os.path.isdir(f):
            force_symlinks_recursive(f)


if len(sys.argv) != 2:
    print("Run: " + sys.argv[0] + " <dir> ")
    exit(1)

force_symlinks_recursive(sys.argv[1])

            
