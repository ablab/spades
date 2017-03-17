#!/usr/bin/env python

import pandas
from pandas import DataFrame

from math import isnan
import os.path
import sys

quast_dir = sys.argv[1]
table_file = sys.argv[2]
heatmap_file = None
if len(sys.argv) > 3:
    heatmap_file = sys.argv[3]

res_table = DataFrame(columns=["bin", "ref", "GF", "purity", "NGA50", "misassemblies"])
gf_table = pandas.read_table(os.path.join(quast_dir, "summary", "TSV", "Genome_fraction_(%).tsv"), dtype=str).set_index("Assemblies")
gfs = gf_table.apply(pandas.to_numeric, errors="coerce")
best_ref = gfs.apply(lambda col: col.idxmax())

for bin, ref in best_ref.iteritems():
    if type(ref) is float:
        row = {"bin": bin, "GF": "-", "ref": "unknown", "purity": "-", "NGA50": "-", "misassemblies": "-"}
    else:
        all_stats = pandas.read_table(os.path.join(quast_dir, "runs_per_reference", ref, "report.tsv"), index_col=0)
        col = all_stats.get(bin)
        purity = 100 - float(col["Unaligned length"]) / float(col["Total length"]) * 100
        row = {"bin": bin, "GF": col["Genome fraction (%)"], "ref": ref, "purity": "{0:.2f}".format(purity),
               "NGA50": col["NGA50"], "misassemblies": col["# misassemblies"]}
    res_table = res_table.append(row, ignore_index=True)

with open(table_file, "w") as out_file:
    res_table.to_csv(out_file, index=False, sep="\t")

if heatmap_file:
    import matplotlib
    # Force matplotlib to not use any Xwindows backend.
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    #gfs = gfs.pivot("reference", "cluster", "GF (%)")
    plot = sns.heatmap(gfs, square=True)
    fig = plot.get_figure()
    fig.savefig(heatmap_file, bbox_inches="tight")
