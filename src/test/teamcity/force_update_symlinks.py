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

def force_symlinks_for_suffix(contig_dir, suffix, name, branch_name):
    contigs = sorted(glob.glob(os.path.join(contig_dir, "2*" + suffix + "_" + branch_name + ".fasta")))
    if len(contigs) == 0:
        return
    clink = os.path.join(contig_dir, name)
    if os.path.islink(clink):
        os.remove(clink)
    os.symlink(os.path.basename(contigs[-1]), clink)
    print("Forced symlink to " + contigs[-1] + " -> " + clink)


def force_symlinks(contig_dir, branch_name):
    force_symlinks_for_suffix(contig_dir, "contigs", "latest_contigs.fasta", branch_name)
    force_symlinks_for_suffix(contig_dir, "scaffolds", "latest_scaffolds.fasta", branch_name)
    force_symlinks_for_suffix(contig_dir, "before_rr", "latest_before_rr.fasta", branch_name)
    force_symlinks_for_suffix(contig_dir, "preliminary", "latest_preliminary.fasta", branch_name)
    force_symlinks_for_suffix(contig_dir, "transcripts", "latest_transcripts.fasta", branch_name)


def force_symlinks_recursive(storage_dir, branch_name):
    force_symlinks(storage_dir, branch_name)
    files = sorted(glob.glob(os.path.join(storage_dir, "*")))
    for f in files:
        if os.path.isdir(f):
            force_symlinks_recursive(f, branch_name)


if len(sys.argv) < 2 or len(sys.argv) > 3:
    print("Run: " + sys.argv[0] + " <dir> [branch name]")
    exit(1)

branch_name = "master" if len(sys.argv) == 2 else sys.argv[2]
force_symlinks_recursive(sys.argv[1], branch_name)

            
