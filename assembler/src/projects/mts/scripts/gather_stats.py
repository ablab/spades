#!/usr/bin/env python
from __future__ import (print_function)

import argparse
from collections import OrderedDict
import os.path
from operator import add
import shutil
import sys

import matplotlib.colors
from sklearn.cluster.bicluster import SpectralCoclustering
import numpy as np
import pandas
from pandas import DataFrame
import json

def bicluster(*cotables):
    table = cotables[0]
    model = SpectralCoclustering(n_clusters=table.shape[1], random_state=0)
    model.fit(table.as_matrix())
    return [cotable.iloc[np.argsort(model.row_labels_), np.argsort(model.column_labels_)] for cotable in cotables]

parser = argparse.ArgumentParser(description="MTS stats extractor")
parser.add_argument("dir", type=str, help="QUAST output directory")
parser.add_argument("name", type=str, help="Output base name")
parser.add_argument("--gf", action="store_true", help="Genome fraction")
parser.add_argument("--nga", action="store_true", help="NGA50 plots for references")
parser.add_argument("--problematic", action="store_true", help="Problematic references report")
parser.add_argument("--summary", action="store_true", help="Produce summary table")
parser.add_argument("--plot", action="store_true", help="Draw plots for metrics")

args = parser.parse_args()

# Load genome fraction table
def table_path(name, combined=False):
    if combined:
        return os.path.join(args.dir, "combined_reference", name)
    return os.path.join(args.dir, "summary", "TSV", name)

gf_table = pandas.read_table(table_path("Genome_fraction_(%).tsv"), index_col=0, na_values="-")
gf_table.index = gf_table.index.map(str)
gf_table.index.name = "ref"
gf_table.fillna(0, inplace=True)

#Drop _broken scaffolds from MetaQUAST report
gf_table.drop([col for col in gf_table.columns if col.endswith("_broken")], axis=1, inplace=True)

# Drop bad bins and missing references
EPS = 0.01
presented_refs = gf_table.apply(lambda row: row.sum() > EPS, axis=1)
lost_refs = gf_table.index[~presented_refs]
gf_table = gf_table.loc[presented_refs, gf_table.apply(lambda col: col.sum() > EPS)]

# Load purity table
total_lengths = pandas.to_numeric(pandas.read_table(table_path("report.tsv", True), index_col=0).loc["Total length"])

aligned_table = pandas.read_table(table_path("Total_aligned_length.tsv"), index_col=0, na_values="-")
aligned_table.fillna(0, inplace=True)

purity_table = aligned_table.apply(lambda col: col / total_lengths[col.name] * 100.)
purity_table.index = purity_table.index.map(str)
purity_table = purity_table[gf_table.columns]

# Detect correspondence between bins and their best references

#gfs = gf_table.loc[presented_refs, gfs.apply(lambda col: col.sum() > EPS)]

best_ref = gf_table.apply(lambda col: col.idxmax())
best_bin = gf_table.apply(lambda row: row.idxmax(), axis=1)
with open(args.name + "_best.tsv", "w") as out_file:
    best_ref.to_csv(out_file, sep="\t")

if args.summary:
    # Load misassemblies from combined reference
    with open(table_path("report.html", True)) as input:
        read = False
        for line in input:
            if "<div id='total-report-json'>" in line:
                read = True
                continue
            if read:
                report = json.loads(line)
                break
        data = [sub[1][1][0]["values"] for sub in report["subreports"]]
        mis_table = pandas.DataFrame(data)
        mis_table.index = report["subreferences"]
        mis_table.columns = report["assembliesNames"]

    # Prepare the summary table
    res_dict = OrderedDict()
    for ref, bin in best_bin.iteritems():
        if type(bin) is float:
            row = ["unknown", "-", "-", "-", "-"]
        else:
            all_stats = pandas.read_table(os.path.join(args.dir, "runs_per_reference", ref, "report.tsv"), index_col=0)
            col = all_stats.get(bin)
            if col is None:
                print("WRONG:", bin, ref)
            purity = purity_table.get_value(ref, bin)
            row = [bin, col["Genome fraction (%)"], purity, col["NGA50"], mis_table.loc[ref, bin]]
        res_dict[ref] = row
    res_table = pandas.DataFrame.from_dict(res_dict, "index")
    res_table.columns = ["bin", "GF", "purity", "NGA50", "misassemblies"]
    res_table.index.name = "ref"
    res_table.to_csv(args.name + "_summary.tsv", index=True, sep="\t", float_format="%.2f")

if args.plot:
    try:
        import matplotlib
        # Force matplotlib to not use any Xwindows backend.
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sns

        # Column permutation to diagonalize the table saving the reference order
        def best_permutation(table):
            rows = list(table.index)
            unused_cols = set(table.columns)
            res = []
            for i in range(len(rows)):
                row = rows[i]
                all_options = [(col, table.get_value(row, col)) for col in unused_cols]
                all_options.sort(key=lambda x: x[1], reverse=True)
                options = list(filter(lambda v: v[1] > 10, all_options)) # (table.loc[rows[i+1:], col] < 0.1).any()]
                #if not options and all_options:
                #    options = [all_options[0]]
                new = [k for k, _ in options]
                res += new
                unused_cols.difference_update(new)
            return res

        def color(v1, v2):
            return matplotlib.colors.hsv_to_rgb((v2 / 200. + 0.16, v1 / 100., 1))

        def combine_with(t1, t2, f):
            return t1.apply(lambda k: [f(x, y) for x, y in zip(k, t2[k.name])])

        def diag_heatmap(t1, t2, filename, best=None):
            perm = best_permutation(t1)
            t1_aligned = t1[perm]
            t2_aligned = t2.loc[t1_aligned.index, perm]

            colors_aligned = combine_with(t1_aligned, t2_aligned, color)

            fig, ax = plt.subplots()

            #TODO: cell text
            table = plt.table(rowLabels=t1_aligned.index, colLabels=t1_aligned.columns,
                              cellColours=colors_aligned.as_matrix(), loc="center",
                              bbox=(0,0,1,1))
            table.auto_set_font_size(False)
            table.set_fontsize(12)
            (bx, by) = t1_aligned.shape
            def unselect(cell):
                cell.set_linewidth(0)

            def select(cell):
                cell.set_linewidth(2)

            for (x, y), cell in table.get_celld().items():
                unselect(cell)
                x -= 1 #Cells indices in the table seems to be shifted by one row
                if best is None or x < 0 or x >= bx or y < 0 or y >= by:
                    continue
                if best[t1_aligned.index[x]] == t1_aligned.columns[y]:
                    select(cell)

            #table.scale(0.4, 2)
            plt.axis("off")

            plt.savefig(filename, bbox_inches="tight")
            plt.clf()

    except:
        print("Cannot import matplotlib and/or seaborn; drawing will be disabled")
        args.plot = False

# GF
if args.gf:
    gf_table.to_csv(args.name + "_gf.tsv", index=True, sep="\t", float_format="%.2f")

    # Draw GF heatmap
    if args.plot:
        diag_heatmap(gf_table, purity_table, args.name + "_gf.png", best_bin)

# Draw NGA plots
if args.nga:
    nga_table = pandas.read_table(os.path.join(args.dir, "summary/TSV/NGA50.tsv"), index_col=0)
    nga_table = nga_table.apply(pandas.to_numeric, errors="coerce")
    nga_series = nga_table.max(axis=1)
    nga_series.to_csv(args.name + "_nga.tsv", sep="\t")
    if args.plot:
        nga_series.index.name = "Reference"
        plot = nga_series.plot(kind="bar", logy=True)
        plot.set_ylabel("Contig length (bp)")
        fig = plot.get_figure()
        fig.savefig(args.name + "_nga.png", bbox_inches="tight")
        plt.clf()

# (Optional) Write summary for problematic references
if args.problematic:
    BAD_THRESHOLD = 90
    ZERO_THRESHOLD = 5
    total_gf_ref = gf_table.sum(1)
    max_gf_ref = gf_table.max(1)
    nonzeroes = gf_table.applymap(lambda x: x > ZERO_THRESHOLD)
    nonzeroes_cnt_ref = nonzeroes.sum(1)
    good_refs = list()
    with open(args.name + "_problems.txt", "w") as out_file:
        for ref in lost_refs:
            print(ref, "is missing completely", file=out_file)
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
    if args.plot:
        drop_refs = [ref for ref in good_refs if best_bin[ref] in good_bins]
        drop_bins = [bin for bin in good_bins if best_ref[bin] in good_refs]
        bad_table = gf_table.drop(drop_refs, axis=0).drop(drop_bins, axis=1)
        bad_table = bad_table[bad_table.apply(lambda row: row.sum() > EPS, axis=1)]
        if bad_table.size:
            diag_heatmap(bad_table, purity_table, args.name + "_bad.png")
