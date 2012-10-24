#!/usr/bin/env python

import os
import sys

sys.path.append(os.path.join(os.path.abspath(sys.path[0]), '../../../src/spades_pipeline'))
from process_cfg import *

def process(config, dataset):
    subst = dict()
    proj = os.path.basename(dataset)
    proj = os.path.splitext(proj)[0]

    subst["project"] = proj
    ds = load_config_from_file(dataset)
    subst["reference"] = ""
    subst["genes"] = ""
    subst["operons"] = ""
    if "reference_genome" in ds.__dict__:
        subst["reference"] = ds.reference_genome
    if "E.coli" in subst["reference"]:
        subst["genes"] = "./data/input/E.coli/genes/genes.txt"
        subst["operons"] = "./data/input/E.coli/genes/operons.txt"

    gendir = os.path.splitext(config)[0]
    if not os.path.exists(gendir):
        os.makedirs(gendir)
    f = open(os.path.join(gendir, dataset), "w")
    for line in file_lines(config):
        for (key, value) in subst.items():
            line = line.replace("$" + key, value)
        f.write(line)
    f.close()


def main():
    for ds in sys.argv[2:]:
        process(sys.argv[1], ds)

if __name__ == '__main__':
    main()
