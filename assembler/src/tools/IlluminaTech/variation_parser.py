############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import os
cur_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(cur_dir, "..", "..", "spades_pipeline")))
import SeqIO
import re

__author__ = 'anton'

class VariationInfo:
    def __init__(self, chr, start, end, len, allele):
        self.chr = chr
        self.start = start
        self.end = end
        self.len = len
        self.allele = allele

    def InSegment(self, pos):
        return pos >= self.start - 2000 and pos <= self.end + 2000

    def SegmentIntersect(self, start, end):
        return self.InSegment(start) or self.InSegment(end) or (self.start > start and self.end < end)

    def Intersect(self, chr, start, end):
        return chr == self.chr and self.SegmentIntersect(start, end)

class VariationParser:
    def __init__(self, file_name):
        self.file = file_name

    def extract(self, data, pattern):
        for s in re.split("[\";]", data):
            if s.startswith(pattern):
                return int(s[len(pattern):])
        return None

    def parse(self):
        f = open(self.file, "r")
        for data in [s.strip().split("\t") for s in f.readlines()][1:]:
            yield VariationInfo(data[0], int(data[1]), self.extract(data[7], "END="), self.extract(data[7], "SVLEN="), data[10] == "1/1")


class VariationFinder:
    chrarr = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chrX", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr20", "chrY", "chr19", "chr22", "chr21", "chr6_ssto_hap7", "chr6_mcf_hap5", "chr6_cox_hap2", "chr6_mann_hap4", "chr6_apd_hap1", "chr6_qbl_hap6", "chr6_dbb_hap3", "chr17_ctg5_hap1", "chr4_ctg9_hap1", "chr1_gl000192_random", "chrUn_gl000225", "chr4_gl000194_random", "chr4_gl000193_random", "chr9_gl000200_random", "chrUn_gl000222", "chrUn_gl000212", "chr7_gl000195_random", "chrUn_gl000223", "chrUn_gl000224", "chrUn_gl000219", "chr17_gl000205_random", "chrUn_gl000215", "chrUn_gl000216", "chrUn_gl000217", "chr9_gl000199_random", "chrUn_gl000211", "chrUn_gl000213", "chrUn_gl000220", "chrUn_gl000218", "chr19_gl000209_random", "chrUn_gl000221", "chrUn_gl000214", "chrUn_gl000228", "chrUn_gl000227", "chr1_gl000191_random", "chr19_gl000208_random", "chr9_gl000198_random", "chr17_gl000204_random", "chrUn_gl000233", "chrUn_gl000237", "chrUn_gl000230", "chrUn_gl000242", "chrUn_gl000243", "chrUn_gl000241", "chrUn_gl000236", "chrUn_gl000240", "chr17_gl000206_random", "chrUn_gl000232", "chrUn_gl000234", "chr11_gl000202_random", "chrUn_gl000238", "chrUn_gl000244", "chrUn_gl000248", "chr8_gl000196_random", "chrUn_gl000249", "chrUn_gl000246", "chr17_gl000203_random", "chr8_gl000197_random", "chrUn_gl000245", "chrUn_gl000247", "chr9_gl000201_random", "chrUn_gl000235", "chrUn_gl000239", "chr21_gl000210_random", "chrUn_gl000231", "chrUn_gl000229", "chrM", "chrUn_gl000226", "chr18_gl000207_random"]
    def __init__(self, file):
        self.variations = list(VariationParser(file).parse())

    def find(self, reference_dir):
        files = [os.path.join(reference_dir, file) for file in os.listdir(reference_dir) if os.path.isfile(os.path.join(reference_dir, file)) and file.endswith("fasta")]
        for file in files:
            sys.stdout.write("Processing file " + file + "\n")
            for rec in SeqIO.parse_fasta(open(file, "r")):
                lines = filter(None, re.split("[\-()]", rec.id))
#                print lines
                chr = lines[0][3:-1]
                left = int(lines[1])
                right = int(lines[2])
#                print chr, left, right
                for variation in self.variations:
                    if variation.Intersect(chr, left, right):
                        sys.stdout.write("Segment chr_{0}_({1}-{2}) intersected with variation ({3}, {4}, {5}, {6}) \n".format(chr, left, right, variation.chr, variation.start, variation.end, variation.allele))

VariationFinder(sys.argv[1]).find(sys.argv[2])
