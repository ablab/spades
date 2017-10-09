#!/usr/bin/env python
from __future__ import (print_function)

import argparse
import os
import os.path

import numpy as np
import pandas
from pandas import DataFrame, Series

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description="MultiStats plotter")
parser.add_argument("--runs", type=str, nargs="+")
parser.add_argument("--plot", action="store_true", help="Draw comparative plots")
parser.add_argument("dir", type=str, help="Multirun output directory")

args = parser.parse_args()

series = list()
if not args.runs:
    args.runs = [run for run in sorted(os.listdir(args.dir)) if "_" in run]

tables = dict()
order = None
for rundir in args.runs:
    params = rundir.split("_")
    assembler = params[0]
    binner = params[1]
    last = "reassembly" if assembler == "main" else "binning"
    try:
        summary_path = os.path.join(args.dir, rundir, "stats", "summary", last + "_summary.tsv")
        summary_table = pandas.read_table(summary_path, index_col="ref", dtype=str)
        nga_path = os.path.join(args.dir, rundir, "stats", "summary", last + "_nga.tsv")
        nga_series = pandas.read_table(nga_path, names=["ref", "NGA50"], index_col=0, dtype=str, squeeze=True)
    except:
        print("Cannot process stats from", summary_path, "; skipping")
        continue
    summary_table["NGA50"] = nga_series
    if not order:
        order = list(summary_table.columns)
        order.remove("bin")

    tables[rundir] = summary_table.stack()

big_table = pandas.concat(tables, axis=1)
big_table.index.names = ["ref", "metrics"]
big_table = big_table.reindex(order, level=1)
big_table.index.drop_duplicates()
big_table.to_csv("summary/summary.tsv", sep="\t")

for table_name in ["NGA50"]:
    table = big_table.xs(table_name, level=1)
    table = table.apply(pandas.to_numeric, errors="coerce")

    table.sort_values(by=table.columns[0], inplace=True)
    fig, ax = plt.subplots()
    for conf, marker in zip(tables.keys(), ".ovx+*D"):
        table[conf].plot(kind="line", logy=True, rot=45, marker=marker, ax=ax)
    ax.legend()

    fig.savefig("summary/" + table_name + ".png", bbox_inches="tight")
    plt.gcf().clear()
