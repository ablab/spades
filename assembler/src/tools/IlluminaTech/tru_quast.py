############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################
import itertools

__author__ = 'anton'

import sys
import os
sys.path.append("src/spades_pipeline")
sys.path.append("src/spades_pipeline/common")
sys.path.append("src/spades_pipeline/truspades")
import barcode_extraction
import parallel_launcher
import support
    # names = ['index', 'Assembly', '# contigs (>= 0 bp)', '# contigs (>= 1000 bp)', 'Total length (>= 0 bp)', 'Total length (>= 1000 bp)', '# contigs', 'Largest contig', 'Total length', 'Reference length', 'GC (%)', 'Reference GC (%)', 'N50', 'NG50', 'N75', 'NG75', 'L50', 'LG50', 'L75', 'LG75', '# misassemblies', '# misassembled contigs', 'Misassembled contigs length', '# local misassemblies', '# unaligned contigs', 'Unaligned length', 'Genome fraction (%)', 'Duplication ratio', "# N's per 100 kbp", '# mismatches per 100 kbp', '# indels per 100 kbp', 'Largest alignment', 'NA50', 'NGA50', 'NA75', 'LA50', 'LGA50', 'LA75']
    #names = ['# misassemblies', '# misassembled contigs', 'Misassembled contigs length', '# local misassemblies', 'Unaligned length', 'Genome fraction (%)', '# mismatches per 100 kbp', '# indels per 100 kbp']

def DefaultZero(d, key):
    if key in d:
        return d[key]
    else:
        return "0"


def RunBarcodeQuast(barcodes, barcode_quast_dir, reference_dir, threads):
    quast_format = " ".join(["quast", "--min-contig", "1000", "--contig-thresholds", "5000,8000,12000", "-e", "-R",
            os.path.join(reference_dir, "{0}.fasta"),
            "{1}", "-o",
            os.path.join(barcode_quast_dir, "{0}")])
    commands = [(barcode_id, quast_format.format(barcode_id, file_name)) for (barcode_id, file_name) in barcodes]
    support.recreate_dir(barcode_quast_dir)
    task = parallel_launcher.ExternalCallTask("", "")
    errors = parallel_launcher.run_in_parallel(task, commands, threads)
    sys.stderr.write(str(errors) + " barcodes failed")


def ParseResults(barcode_quast_dir, ids):
    names = []
    reports = list()
    for barcode_id in ids:
        report_path = os.path.join(barcode_quast_dir, barcode_id, "report.tsv")
        new_item = dict()
        new_item["id"] = barcode_id
        if os.path.isfile(report_path):
            report = open(report_path, "r")
            for param in report:
                line = param.split("\t")
                new_item[line[0]] = line[1].strip()
                if not line[0] in names:
                    names.append(line[0])
        reports.append(new_item)
    return names, reports


def CollectResults(names, reports):
    values = {metric: 0 for metric in names}
    values["#partially unaligned"] = 0
    for line in reports:
        for name in names:
            if name != "# unaligned contigs" and name != "index" and name != "Assembly":
                values[name] += float(DefaultZero(line, name))
            else:
                tmp = DefaultZero(line, name).split(" ")
                if len(tmp) > 1:
                    values["# unaligned contigs"] += int(tmp[0])
                if len(tmp) >= 3:
                    values["#partially unaligned"] += int(tmp[2])
    return values


def RunTruQuast(input_dir, reference_dir, output_dir, threads):
    support.ensure_dir_existence(output_dir)
    if os.path.exists(os.path.join(input_dir, "dataset.info")):
        ids = [barcode.id for barcode in barcode_extraction.ReadDataset(os.path.join(input_dir, "dataset.info"))]
        files = [os.path.join(input_dir, "barcodes", bid, "truseq_long_reads.fasta") for bid in ids]
    else:
        files = [os.path.join(input_dir, file) for file in os.listdir(input_dir) if file.endswith(".fasta") or file.endswith(".fa")]
        ids = [f[:f.rfind(".")] for f in os.listdir(input_dir) if file.endswith(".fasta") or file.endswith(".fa")]

    barcode_quast_dir = os.path.join(output_dir, "barcode_quast")
    RunBarcodeQuast(zip(ids, files), barcode_quast_dir, reference_dir, threads)
    names, reports = ParseResults(barcode_quast_dir, ids)
    names.append("#partially unaligned")
    values = CollectResults(names, reports)
    results = open(os.path.join(output_dir, "results.tsv"), "w")
    for name in names:
        results.write(name + "\t" + str(int(values[name])) + "\t" + str(values[name] / len(ids)) + "\n")
    results.close()

if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.stderr.write("Usage: python src/tools/IlluminaTech/tru_quast.py <input_dir> <references_dir> <output_dir> <threads>\n")
    else:
        RunTruQuast(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
