#!/usr/bin/env python
from __future__ import (print_function)

import pandas
from pandas import DataFrame

from math import isnan
import argparse
import os.path
from operator import add
import sys

parser = argparse.ArgumentParser(description="MTS - Metagenomic Time Series")
parser.add_argument("dir", type=str, help="QUAST output directory")
parser.add_argument("name", type=str, help="Output base name")
parser.add_argument("--problematic", action="store_true", help="Problematic references report")
parser.add_argument("--heatmap", action="store_true", help="Best reference summary table")

args = parser.parse_args()

# Write summary table with correspondence between bins and their best references
res_table = DataFrame(columns=["bin", "ref", "GF", "purity", "NGA50", "misassemblies"])
gf_table = pandas.read_table(os.path.join(args.dir, "summary", "TSV", "Genome_fraction_(%).tsv"), dtype=str).set_index("Assemblies")
gfs = gf_table.apply(pandas.to_numeric, errors="coerce")
best_ref = gfs.apply(lambda col: col.idxmax())

for bin, ref in best_ref.iteritems():
    if type(ref) is float:
        row = {"bin": bin, "GF": "-", "ref": "unknown", "purity": "-", "NGA50": "-", "misassemblies": "-"}
    else:
        all_stats = pandas.read_table(os.path.join(args.dir, "runs_per_reference", ref, "report.tsv"), index_col=0)
        col = all_stats.get(bin)
        purity = 100 - float(col["Unaligned length"]) / float(col["Total length"]) * 100
        row = {"bin": bin, "GF": col["Genome fraction (%)"], "ref": ref, "purity": "{0:.2f}".format(purity),
               "NGA50": col["NGA50"], "misassemblies": col["# misassemblies"]}
    res_table = res_table.append(row, ignore_index=True)

with open(args.name + "_summary.tsv", "w") as out_file:
    res_table.to_csv(out_file, index=False, sep="\t")

# (Optional) Draw GF heatmap
if args.heatmap:
    try:
        import matplotlib
        # Force matplotlib to not use any Xwindows backend.
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sns

        plot = sns.heatmap(gfs, square=True)
        fig = plot.get_figure()
        fig.savefig(args.name + "_gf.png", bbox_inches="tight")
        plt.gcf().clear()
    except:
        print("Can't import matplotlib and/or seaborn; heatmap drawing will be disabled")
        args.heatmap = False

# (Optional) Write summary for problematic references
if args.problematic:
    BAD_THRESHOLD = 90
    ZERO_THRESHOLD = 5
    total_gf_ref = gfs.sum(1)
    max_gf_ref = gfs.max(1)
    nonzeroes = gfs.applymap(lambda x: x > ZERO_THRESHOLD)
    nonzeroes_cnt_ref = nonzeroes.sum(1)
    good_refs = list()
    with open(args.name + "_problems.txt", "w") as out_file:
        for ref, gf in total_gf_ref.iteritems():
            if max_gf_ref[ref] < BAD_THRESHOLD:
                if gf < BAD_THRESHOLD:
                    print(ref, "is underassembled: at least", 100 - gf, "% GF was lost", file=out_file)
                else:
                    print(ref, "is fractured: best bin is only", max_gf_ref[ref], "% GF", file=out_file)
                continue
            if nonzeroes_cnt_ref[ref] > 1:
                print(ref, "is presented in", nonzeroes_cnt_ref[ref], "bins", file=out_file)
                continue
            good_refs.append(ref)
        nonzeroes_cnt_bin = nonzeroes.sum(0)
        good_bins = list()
        for bin, cnt in nonzeroes_cnt_bin.iteritems():
            if cnt > 1:
                print(bin, "is a mixture of", cnt, "references", file=out_file) #TODO: which ones?
            else:
                good_bins.append(bin)
    if args.heatmap:
        bad_table = gfs.drop(good_refs, axis=0).drop(good_bins, axis=1)
        if bad_table.size:
            plot = sns.heatmap(bad_table, square=True)
            fig = plot.get_figure()
            fig.savefig(args.name + "_bad.png", bbox_inches="tight")
