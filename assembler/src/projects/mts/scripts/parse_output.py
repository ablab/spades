#!/usr/bin/env python
from __future__ import print_function

import argparse
import os.path

argparser = argparse.ArgumentParser(description="Binner output formatter")
argparser.add_argument("--type", "-t", type=str, help="Binner type (canopy or concoct)", default="canopy")
argparser.add_argument("--output", "-o", type=str, help="Output directory with annotations")
argparser.add_argument("input", type=str, help="File with binning info")

class Parser:
    def __init__(self):
        self.samples_annotation = dict()

    def add(self, line):
        sample_contig, bin_id = self.parse(line)
        sample_contig = sample_contig.split('-', 1)
        sample = sample_contig[0]
        contig = sample_contig[1]
        if sample not in self.samples_annotation:
            self.samples_annotation[sample] = dict()

        annotation = self.samples_annotation[sample]
        if contig not in annotation:
            annotation[contig] = list()

        annotation[contig].append(bin_id)

class CanopyParser(Parser):
    def parse(self, line):
        annotation_str = line.split()
        bin_id = annotation_str[0].strip()
        sample_contig = annotation_str[1].strip()
        return (sample_contig, bin_id)

class ConcoctParser(Parser):
    def parse(self, line):
        annotation_str = line.split(",", 1)
        bin_id = annotation_str[1].strip()
        sample_contig = annotation_str[0].replace("~", ",")
        return (sample_contig, bin_id)

parsers = {"canopy": CanopyParser(), "concoct": ConcoctParser()}

args = argparser.parse_args()
parser = parsers[args.type]

with open(args.input, "r") as input_file:
    for line in input_file:
        parser.add(line)

for sample, annotation in parser.samples_annotation.items():
    with open(os.path.join(args.output, sample + ".ann"), "w") as sample_out:
        annotation = parser.samples_annotation[sample]

        for contig in annotation:
            print(contig, ":", " ".join(annotation[contig]), file=sample_out)
