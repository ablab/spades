#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sam_parser

import SeqIO

import sys

class PatternContigFilter:
    def __init__(self, contigs, sam, pattern, rc_pattern):
        self.sam = sam
        reads = []
        self.filter = [False] * len(contigs)
        for rec in sam:
            reads.append(rec)
            if len(reads) == 2:
                left_sequence = str(reads[0].seq.upper())
                right_sequence = str(reads[1].seq.upper())
                if left_sequence.find(pattern) != -1 or right_sequence.find(rc_pattern) != -1 or right_sequence.find(pattern) != -1 or left_sequence.find(rc_pattern) != -1:
                    if not reads[0].is_unmapped:
                        self.filter[reads[0].tid] = True
                    if not reads[1].is_unmapped:
                        self.filter[reads[1].tid] = True
                reads = []

    def Filter(self, contig):
        return self.filter[self.sam.gettid(contig.id)]

class ContigLengthFilter:
    def __init__(self, min_length):
        self.min_length = min_length

    def Filter(self, contig):
        return len(contig.seq) >= self.min_length

# def Filter(contigs, left_reads, right_reads, length_threshold, pattern):
#     dataset_data = pyyaml.load(yaml_file, 'r')
