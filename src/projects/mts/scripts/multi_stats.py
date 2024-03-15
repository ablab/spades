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
matplotlib.rcParams["axes.labelsize"] = "small"

import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description="MultiStats plotter")
parser.add_argument("--runs", type=str, nargs="+")
parser.add_argument("--assemblers", type=str, nargs="+", default=["main", "groups", "megahit", "spades"])
parser.add_argument("--binners", type=str, nargs="+", default=["concoct"])
parser.add_argument("--plot", action="store_true", help="Draw comparative plots")
parser.add_argument("-r", "--rename", action="store_true")
parser.add_argument("--strain", type=str, default=None, help="Keep a single strain only")
parser.add_argument("--dir", type=str, default=".", help="Multirun output directory")

args = parser.parse_args()
args.output = "summary"

if not os.path.isdir(args.output):
    os.mkdir(args.output)

def fullname(name):
    return os.path.join(args.output, name)

series = list()
if not args.runs:
    args.runs = ["{}_{}".format(assembler, binner) for assembler in args.assemblers for binner in args.binners]
    if not args.runs:
        args.runs = [run for run in sorted(os.listdir(args.dir)) if run.count("_") == 1]

def rename_ref(name):
    name = name.replace("-", "_")
    params = name.split("_")
    if len(params) > 1:
        params[1] = params[1][:3]
    else:
        params[0] = params[0][:5]
    return "_".join(params)

tables = dict()
order = None
table_names = list()
for rundir in args.runs:
    params = rundir.split("_")
    assembler = params[0]
    binner = params[1]
    last = "reassembly" if assembler == "main" or assembler == "groups" else "binning"
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
    table_name = rundir
    if args.rename:
        rename_dict = {"main": "mts_default", "groups": "mts_groups", "spades": "alt_metaspades", "megahit": "alt_megahit"}
        table_name = rename_dict[assembler]
    summary_table.index = map(rename_ref, summary_table.index)
    if args.strain:
        def skip_strain(ref):
            return ref.count("_") > 1 and not ref.endswith(args.strain)
        summary_table.drop([ref for ref in summary_table.index if skip_strain(ref)], axis=0, inplace=True)
    summary_table.sort_index(inplace=True)
    tables[table_name] = summary_table.stack()
    table_names.append(table_name)

big_table = pandas.concat(tables, axis=1, keys=table_names)
big_table.index.names = ["ref", "metrics"]
big_table = big_table.reindex(order, level=1)
big_table.index.drop_duplicates()
big_table.to_csv(fullname("summary.tsv"), sep="\t")

stat_labels = {"GF": "GF (%)", "purity": "purity (%)", "NGA50": "NGA50 (bp)", "misassemblies": "misasm"}
stat_pos = {"GF": (0, 0), "purity": (0, 1), "NGA50": (1, 0), "misassemblies": (1, 1)}

fig, combo = plt.subplots(2, 2, sharex=True)

for table_name in stat_labels:
    table = big_table.xs(table_name, level=1)
    table = table.apply(pandas.to_numeric, errors="coerce")

    i, j = stat_pos[table_name]

    sfig, single = plt.subplots()
    # for conf, marker in zip(tables.keys(), ".ovx+*D"):
    #     table[conf].plot(kind="line", logy=True, rot=45, marker=marker, ax=ax)
    # ax.legend()

    logy = table_name == "NGA50"

    for ax in combo[i, j], single:
        table.plot(ax=ax, kind="bar", logy=logy, rot=90)
        ax.legend_.remove()
        ax.set(xlabel="")

    #plt.ylabel(stat_labels[table_name])

    sfig.savefig(fullname(table_name + ".png"), bbox_inches="tight")

handles, labels = combo[0, 0].get_legend_handles_labels()
fig.legend(handles, labels , loc="upper center", ncol=len(args.runs))
fig.savefig(fullname("combo.png"), bbox_inches="tight", dpi=300)

def color(v1, v2):
    return matplotlib.colors.hsv_to_rgb((v1 / 200. + 0.16, v2 / 100., 1))

def combine_with(t1, t2, f):
    return t1.apply(lambda k: [f(x, y) for x, y in zip(k, t2[k.name])])

def packed_heatmap(t1, t2, filename, best=None):
    t1_aligned = t1
    t2_aligned = t2

    colors_aligned = combine_with(t1_aligned, t2_aligned, color)

    w, h = t1.shape

    font_size = 12
    cell_height = font_size / 72.0 * 3
    cell_width = 0.6
    fig = plt.figure(figsize=(1 + cell_width * w, cell_height * (1 + h)))

    #TODO: cell text
    table = plt.table(rowLabels=t1_aligned.index, colLabels=t1_aligned.columns,
                      cellColours=colors_aligned.as_matrix(), loc="center",
                      bbox=(0,0,1,1))
    table.auto_set_font_size(False)
    table.set_fontsize(12)

    #table.scale(0.4, 2)
    plt.axis("off")

    plt.savefig(filename, bbox_inches="tight")
    plt.clf()

gf_table = big_table.xs("GF", level=1).apply(pandas.to_numeric)
purity_table = big_table.xs("purity", level=1).apply(pandas.to_numeric)
packed_heatmap(gf_table, purity_table, fullname("GF_purity.png"))
