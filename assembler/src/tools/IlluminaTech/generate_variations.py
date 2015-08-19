############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import random
from Bio.Graphics.BasicChromosome import _ChromosomeComponent
import itertools

__author__ = 'anton'
import sys

nucls = "ACGT"
default_rate = 0.05
len_range = (3000, 15000)
sys.path.append("src/spades_pipeline")

import SeqIO


def Mutate(seq, rate):
    s = list(seq)
    for i in range(len(seq)):
        if s[i] in nucls and random.random() < rate:
            s[i] = random.choice(nucls)
    return "".join(s)


def RandSegment(reference, len_r):
    pos = random.randint(0, len(reference) - len_r[1])
    return reference[pos : pos + random.randint(len_r[0], len_r[1])]


def GenerateInsertions(numins, result):
    ref = "".join(result)
    insertions = []
    for i in range(numins):
        seq = Mutate(RandSegment(ref, len_range), default_rate)
        insertions.append((random.randint(0, len(result)), seq))
    return sorted(insertions)


def GenerateDeletions(numdel, result):
    deletions = []
    for i in range(numdel):
        l = random.randint(0, len(result) - len_range[1])
        r = l + random.randint(len_range[0], len_range[1])
        if result[l:r].find("$") == -1:
            deletions.append((l, r - l))
    return sorted(deletions)

def GroupByChrom(positions, reference):
    last = 0
    i = 0
    result = []
    for l in [len(r) for r in reference]:
        last += l
        tmp = []
        while i < len(positions) and positions[i][0] < last:
            tmp.append((positions[i][0] + l - last, positions[i][1]))
            i += 1
        result.append(tmp)
    return result

def Apply(seq, ins, d):
    result = []
    i = 0
    j = 0
    last = 0
    l = 0
    while i < len(ins) or j < len(d):
        if i < len(ins) and (j == len(d) or ins[i][0] < d[j][0]):
            if last < ins[i][0]:
                result.append(seq[last:ins[i][0]])
                l += ins[i][0] - last
                sys.stdout.write("Insertion: " + str(l) + " " + str(l + len(ins[i][1])) + "\n")
                result.append(ins[i][1])
                l += len(ins[i][1])
                last = ins[i][0]
            i += 1
        else:
            if last < d[j][0]:
                result.append(seq[last:d[j][0]])
                l += d[j][0] - last
                sys.stdout.write("Deletion: " + str(l) + " " + str(d[j][1]) + "\n")
                last = d[j][0] + d[j][1]
            j += 1
    result.append(seq[last:])
    return "".join(result)

def Generate(input, output, numins, numdel):
    reference = list(input)
    result = "".join([ch.seq for ch in reference])
    l = sum([len(ch) for ch in reference])
    ins = GroupByChrom(GenerateInsertions(numins, result), reference)
    d = GroupByChrom(GenerateDeletions(numdel, result), reference)
    for ch_ins, ch_d, chrom in itertools.izip(ins, d, reference):
        sys.stdout.write("Chromosome " + chrom.id + "\n")
        rec = SeqIO.SeqRecord(Apply(chrom.seq, ch_ins, ch_d), chrom.id)
        SeqIO.write(rec, output, "fasta")

if __name__ == '__main__':
    Generate(SeqIO.parse(open(sys.argv[1], "r"), "fasta"), open(sys.argv[2], "w"), int(sys.argv[3]), int(sys.argv[3]))
