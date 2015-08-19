#! /usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#input_folder, output_folder
import sys
import os
import getopt
import math

def DefaultZero(d, key):
    if key in d:
        return d[key]
    else:
        return "0"
    
cnt = 0
reports = list()
input_folder = sys.argv[1]
output_folder = sys.argv[2]
index = [file[:-4] for file in os.listdir(input_folder) if file.endswith(".tsv")]
#names = ['index', 'Assembly', '# contigs (>= 0 bp)', '# contigs (>= 1000 bp)', 'Total length (>= 0 bp)', 'Total length (>= 1000 bp)', '# contigs', 'Largest contig', 'Total length', 'Reference length', 'GC (%)', 'Reference GC (%)', 'N50', 'NG50', 'N75', 'NG75', 'L50', 'LG50', 'L75', 'LG75', '# misassemblies', '# misassembled contigs', 'Misassembled contigs length', '# local misassemblies', '# unaligned contigs', 'Unaligned length', 'Genome fraction (%)', 'Duplication ratio', "# N's per 100 kbp", '# mismatches per 100 kbp', '# indels per 100 kbp', 'Largest alignment', 'NA50', 'NGA50', 'NA75', 'LA50', 'LGA50', 'LA75']
#names = ['# misassemblies', '# misassembled contigs', 'Misassembled contigs length', '# local misassemblies', 'Unaligned length', 'Genome fraction (%)', '# mismatches per 100 kbp', '# indels per 100 kbp']
names = []
values = dict()
for name in names:
    values[name] = 0
for l in index:
    cnt += 1
    report_path = os.path.join(input_folder,  + "/" + l.strip() + ".tsv")
    new_item = dict()
    new_item["index"] = l

    if os.path.isfile(report_path):
        report = open(report_path)
        for param in report:
            params = param.split("\t")
            new_item[params[0]] = params[1].strip()
            if not params[0] in names:
                values[params[0]] = 0
                names.append(params[0])
    reports.append(new_item)
index.close()

all_quast = open(os.path.join(output_folder, "table.tsv"), "w")
all_quast.write("\t".join(names))
names.append("#partially unaligned")
values["#partially unaligned"] = 0
for line in reports:
    all_quast.write("\t".join([DefaultZero(line, x) for x in names]))
    all_quast.write("\n")
    for name in names:
        if name != "# unaligned contigs" and name != "index" and name != "Assembly":
            values[name] += float(DefaultZero(line, name))
        else:
            tmp = DefaultZero(line, name).split(" ")
            if len(tmp) > 1:
                values["# unaligned contigs"] += int(tmp[0])
            if len(tmp) >= 3:
                values["#partially unaligned"] += int (tmp[2])


all_quast.close()
results = open(os.path.join(output_folder, "results.tsv"), "w")
for name in names:
    results.write(name + "\t" + str(int(values[name])) + "\t" + str(values[name] / cnt) + "\n")
results.close()
