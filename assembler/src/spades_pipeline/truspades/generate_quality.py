############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import re
import sys
import itertools
import sam_parser

pattern = re.compile('([0-9]*)([MIDNSHP])')

def parse(cigar, len, pos = 0):
    if cigar == "=" :
        for i in range(len):
            yield (i, i + pos)
        return
    if cigar == "X":
        return
    cur = 0
    curr = pos
    for n, c in pattern.findall(cigar):
        if n:
            n = int(n)
        else:
            n = 1
        if c == 'M':
            for i in range(n):
                yield (cur, curr)
                cur += 1
                curr += 1
        elif c == 'DPN':
            curr += n
        elif c in "IS":
            cur += n

def CollectQuality(contigs, sam):
    qual = [[[0,0] for i in range(len(contig))] for contig in contigs]
    for rec in sam:
        if rec.proper_alignment:
            for seq_pos, contig_pos in parse(rec.cigar, rec.alen, rec.pos - 1):
                if rec.seq[seq_pos] == contigs[rec.tid].seq[contig_pos]:
                    qual[rec.tid][contig_pos][1] += 1
                    qual[rec.tid][contig_pos][0] += ord(rec.qual[seq_pos])
    return qual

def CountContigQuality(contigs, qual):
    for i in range(len(contigs)):
        cnt = 0
        qual_list = [chr(33)] * len(contigs[i])
        for pos in range(len(contigs[i])):
            q = qual[i][pos]
            if q[1] != 0:
                qual_list[pos] = chr(q[0] / q[1])
            else:
                cnt += 1
        contigs[i].qual = "".join(qual_list)


def GenerateQuality(contigs, sam):
    qual = CollectQuality(contigs, sam)
    CountContigQuality(contigs, qual)

