import subprocess
import datetime
import os
import sys


if len(sys.argv) < 6:
    print "Usage: prepare_cfg.py <yaml-template> <graph> <reads> <K> <type>"
    exit(-1)

yaml_template = sys.argv[1]
reads = sys.argv[3]
graph = sys.argv[2]
K = sys.argv[4]
tp = sys.argv[5]

def prepare_cfg(yaml_template, params):
    to_print = ""
    with open(yaml_template, "r") as fin:
        for ln in fin.readlines():
            print_ln = ln
            for p in params.keys():
                if ln.startswith(p):
                    print_ln = p + ": " + params[p] + "\n"
                    break
            to_print += print_ln
    print to_print

params = {"path_to_graphfile": graph, "path_to_sequences": reads, "k": K, "data_type": tp}
prepare_cfg(yaml_template, params)

