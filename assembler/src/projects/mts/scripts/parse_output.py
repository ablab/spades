#!/usr/bin/env python
from __future__ import print_function

import argparse
import os.path

argparser = argparse.ArgumentParser(description="Binner output formatter")
argparser.add_argument("--output", "-o", type=str, help="Output directory with annotations")
argparser.add_argument("input", type=str, help="File with binning info")

class Parser:
    def __init__(self):
        self.samples_annotation = dict()

    def add(self, line):
        sample_contig, bin_id = self.parse(line)
        sample_contig = sample_contig.split('-', 1)
        if len(sample_contig) > 1:
            sample = sample_contig[0]
            contig = sample_contig[1]
        else:
            sample = "group1"
            contig = sample_contig[0]
        if sample not in self.samples_annotation:
            self.samples_annotation[sample] = dict()

        annotation = self.samples_annotation[sample]
        if contig not in annotation:
            annotation[contig] = list()

        annotation[contig].append(bin_id)

    def parse_file(self, file):
        with open(file, "r") as input_file:
            for line in input_file:
                self.add(line)

    def parse(self, line):
        annotation_str = line.split("\t", 1)
        bin_id = annotation_str[1].strip()
        sample_contig = annotation_str[0]
        return (sample_contig, bin_id)

if __name__ == "__main__":
    args = argparser.parse_args()
    parser = Parser()

    parser.parse_file(args.input)

    for sample, annotation in parser.samples_annotation.items():
        with open(os.path.join(args.output, sample + ".ann"), "w") as sample_out:
            annotation = parser.samples_annotation[sample]
            for contig in annotation:
                print(contig, "\t", " ".join(annotation[contig]), sep="", file=sample_out)
