#! /usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import SeqIO
from SeqIO import SeqRecord

import sys
import os
import shutil
import sam_parser
import itertools


def ConstructCoverage(sam, contigs, k):
    cov = dict()
    for contig in range(len(contigs)):
        cov[contig] = [0] * (len(contigs[contig]) + 1)
    reads = []
    for rec in sam:
        reads.append(rec)
        if len(reads) == 2:
            if reads[0].proper_alignment:
                if reads[0].pos + k - 1 < reads[1].pos + reads[1].alen - k:
                    cov[reads[0].tid][reads[0].pos + k - 1] += 1
                    cov[reads[0].tid][reads[1].pos + reads[1].alen - k] -= 1
                else:
                    if reads[1].pos + k - 1 < reads[0].pos + reads[0].alen - k:
                        cov[reads[1].tid][reads[1].pos + k - 1] += 1
                        cov[reads[0].tid][reads[0].pos + reads[0].alen - k] -= 1
            reads = []
    return cov

def ConstructCoverageSingle(sam, contigs, k):
    cov = dict()
    for contig in range(len(contigs)):
        cov[contig] = [0] * (len(contigs[contig]) + 1)
    for rec in sam:
        if rec.proper_alignment:
            if rec.pos + k - 1 < rec.pos + rec.alen - k:
                cov[rec.tid][rec.pos + k - 1] += 1
                cov[rec.tid][rec.pos + rec.alen - k] -= 1
    return cov

def OutputHist(cov, contigs, folder):
    if os.path.exists(folder):
        shutil.rmtree(folder)
    os.makedirs(folder)
    for contig in range(len(contigs)):
        f = open(folder + "/" + contigs[contig].id, "w")
        cur = 0
        for i in range(len(cov[contig])):
            cur += cov[contig][i]
            f.write(str(i) + " " + str(cur) + "\n")
        f.close()

def ConstructSimpleCoverage(sam, contigs, k):
    simple_cov = dict()
    for contig in range(len(contigs)):
        simple_cov[contig] = [0] * (len(contigs[contig]) + 1)
    for rec in sam:
        if not rec.is_unmapped:
            simple_cov[rec.tid][rec.pos] += 1
            simple_cov[rec.tid][rec.pos + rec.alen] -= 1
    return simple_cov

def BreakContig(cov, k, min0):
    l = len(cov) - 1
    if l < 2 * k:
        return []
    result = []
    cur = 0
    cur_len0 = 0
    prev_break = 0
    for i in range(l):
        cur += cov[i]
        if cur == 0:
            cur_len0 += 1
        else:
            if cur_len0 == i:
                prev_break = max(0, i - k)
            elif cur_len0 > min0:
                result.append([prev_break, i - cur_len0])
                prev_break = i
            cur_len0 = 0
    result.append([prev_break, min(l, l - cur_len0 + k)])
    return result

class ContigBreaker:
    def __init__(self, contigs, sam, k, min0):
        self.part_list_ = []
        self.contigs = contigs
        self.sam = sam
        cov = ConstructCoverage(self.sam, contigs, k)
        # OutputHist(cov, contigs, "tmp")
        # simple_cov = ConstructSimpleCoverage(sam, k)
        for contig in range(len(contigs)):
            parts = BreakContig(cov[contig], k, min0)
            self.part_list_.append(parts)

    def Break(self, contig):
        result = []
        #print contig.id
        #print self.sam.gettid(contig.id)
        for part in self.part_list_[self.sam.gettid(contig.id)]:
            result.append(contig.subseq(part[0], part[1]))
        return result

    def OutputBroken(self, output_file):
        output = open(output_file, "w")
        for contig in self.contigs:
            for subcontig in self.Break(contig):
                SeqIO.write(subcontig, output, "fasta")
        output.close()

class PatternBreaker:
    def __init__(self, pattern, rc_pattern, max_cut):
        self.pattern = pattern
        self.rc_pattern = rc_pattern
        self.max_cut = max_cut

    def FindLeftPos(self, seq):
        l1 = seq.find(self.pattern)
        l2 = seq.find(self.rc_pattern)
        if l1 == -1:
            l1 = len(seq)
        if l2 == -1:
            l2 = len(seq)
        l = min(l1, l2) + len(self.pattern)
        if l < self.max_cut:
            return l
        else:
            return 0

    def FindRightPos(self, seq):
        l1 = seq.rfind(self.pattern)
        l2 = seq.rfind(self.rc_pattern)
        if l1 == -1:
            l1 = 0
        if l2 == -1:
            l2 = 0
        l = max(l1, l2)
        if l > len(seq) - self.max_cut:
            return l
        else:
            return len(seq)

    def Break(self, contig):
        if len(contig) < 2 * self.max_cut:
            return []
        l,r = self.FindLeftPos(contig.seq), self.FindRightPos(contig.seq)
        return [contig.subseq(l, r)]

class NBreaker:
    def __init__(self, min_N):
        self.min_N = min_N

    def Break(self, contig):
        result = []
        last_break = 0;
        pos = 0
        while(pos < len(contig) and contig[pos] == 'N'):
            pos += 1
        while pos <len(contig):
            rpos = pos
            while(rpos < len(contig) and contig[rpos] == 'N'):
                rpos += 1
            if rpos - pos >= self.min_N:
                result.append(contig.subseq(last_break, pos))
                last_break = rpos
            pos = max(rpos, pos + 1)
        if last_break != len(contig):
            result.append(contig.subseq(last_break, len(contig)))
        return result

#if __name__ == '__main__':
#    ContigBreaker(sys.argv[1], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])).OutputBroken(sys.argv[2])
