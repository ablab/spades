#!/usr/bin/env python
from __future__ import print_function

import argparse
import os.path

argparser = argparse.ArgumentParser(description="Binner output formatter")
argparser.add_argument("--type", "-t", choices=["canopy", "concoct", "maxbin"], help="Binner type", default="canopy")
argparser.add_argument("--output", "-o", type=str, help="Output directory with unified binning results")
argparser.add_argument("input", type=str, help="File with binning info")

class Parser:
    def __init__(self):
        self.bins = []

    def add(self, line):
        sample_contig, bin_id = self.parse(line)
        self.bins.append(sample_contig + "\t" + bin_id)

    def parse_file(self, file):
        with open(file, "r") as input_file:
            for line in input_file:
                self.add(line)

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

parsers = {"canopy": CanopyParser(), "concoct": ConcoctParser(), "maxbin": ConcoctParser()}

if __name__ == "__main__":
    args = argparser.parse_args()
    parser = parsers[args.type]

    parser.parse_file(args.input)

    with open(args.output, "w") as sample_out:
        for sample in parser.bins:
            print(sample, file=sample_out)
