#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import subprocess
import os
import argparse
import os.path
import sys

import pandas as pd

def exit_with_error(message):
    sys.stderr.write("Error: {__message}\n".format(__message=message))
    sys.exit(1)

class Comparator(object):

    @staticmethod
    def improved_count(corrected):
        return corrected.query('levDistance < levDistance_baseline').shape[0]

    @staticmethod
    def comparable_improved_count(corrected):
        return corrected.query('levDistance < levDistance_baseline and comparable==True').shape[0]

    @staticmethod
    def fill_comparable_aligments_flag(corrected):
        corrected.eval('comparable = ((abs(refStart - refStart_baseline) < 15) and (abs(refEnd - refEnd_baseline) < 15))', inplace=True)
        return corrected


    @staticmethod
    def corrupted_count(corrected):
        corrupted = corrected.query('levDistance > levDistance_baseline and comparable==True')
        print("Corrupted")
        print(corrupted.index)
        return corrupted.shape[0]

    @staticmethod
    def full_fixed_count(corrected):
        return corrected.query("hammingDistance == 0 and hammingDistance_baseline != 0").shape[0]

    @staticmethod
    def comparable_full_fixed_count(corrected):
        return corrected.query("hammingDistance == 0 and hammingDistance_baseline != 0 and comparable==True").shape[0]

    @staticmethod
    def stats_filename(path):
        return path + ".stats.tsv.gz"


    def print_count(self,  count):
        return "{__count:.0f} ({__percent:.2f}%)".format(__count=count, __percent=100.0 * count / self.baseline.shape[0])

    def calc_stats(self, path):
        corrected = pd.read_csv(Comparator.stats_filename(path), sep="\t")
        #wtf, ids after correction contains some trash values in suffix
        corrected["id"]=corrected["id"].apply(lambda  x : x.split("_")[0])
        corrected.set_index("id", inplace=True)
        # sum_stats = corrected.sum(axis=0)
        corrected = self.join_with_baseline(corrected)
        corrected = self.fill_comparable_aligments_flag(corrected)
        corrected.eval('gain = (hammingDistance_baseline - hammingDistance) / (hammingDistance_baseline + 1)')
        mean_stats = corrected.query("comparable==True")[["hammingDistance", "levDistance", "gain"]].mean(axis=0)
        print(corrected.shape)
        corrected_dist_one = corrected.query("hammingDistance_baseline == 1")
        print("Uncomparable reads")
        uncomparable = corrected.query('comparable == False')
        print(uncomparable.index)
        comparable_count = corrected.query('comparable == True').shape[0]
        corrupted = self.corrupted_count(corrected)
        return pd.Series({"path" : path,
                          "comparable_count" : self.print_count(comparable_count),
                          "realigned_count" : self.print_count(uncomparable.shape[0]),
                          # "meanHammingDistance" : mean_stats["hammingDistance"],
                          "mean_lev_distance" : "{__dist:.4f} (x{__percent:.2f})".format(__dist=mean_stats["levDistance"], __percent=self.__mean_lev_dist / mean_stats["levDistance"]),
                          # "gain" : mean_stats["gain"],
                          # "total_insertions" : sum_stats["insertions"],
                          # "total_deletions" : sum_stats["deletions"],
                          # "total_mismatch" : sum_stats["mismatch"],
                          # "improved_count" : self.improved_count(corrected),
                          "improved_count" : self.print_count(self.comparable_improved_count(corrected)),
                          "corrupted_count" : self.print_count(corrupted),
                          "realigned_corrupted_count" : self.print_count(corrupted + uncomparable.shape[0]),
                          # "full_fixed_count" : self.full_fixed_count(corrected),
                          "full_fixed_count" : self.print_count(self.comparable_full_fixed_count(corrected)),
                          "one_error_corrupted_count" : self.print_count(self.corrupted_count(corrected_dist_one)),
                          "one_error_full_fixed_count" : self.print_count(self.comparable_full_fixed_count(corrected_dist_one))})

    def join_with_baseline(self, corrected):
        return corrected.join(self.baseline, rsuffix="_baseline")


    def run_calc_stats_task(self, reads_path):
        stats_path = self.stats_filename(reads_path)
        if os.path.isfile(stats_path) and not self.force_recalc:
            return
        cmd = "java -Xmx64G -jar {__comparator_jar} {__reference} {__reads} {__stats_file}".format(__comparator_jar = self.comparator_jar,
                                                                                           __reference=self.reference,
                                                                                           __reads = reads_path,
                                                                                           __stats_file = stats_path)
        subprocess.call(cmd, shell=True)


    def print_baseline_stats(self):
        print("Baseline distance stats:")
        print("Mean distance")
        print(self.baseline[["hammingDistance", "levDistance"]].mean(axis=0))
        self.__mean_lev_dist = self.baseline[["levDistance"]].mean(axis=0)[0]
        print("Error sums:")
        print(self.baseline.drop(['levDistance', 'hammingDistance'], axis=1).sum())

    def __init__(self, reference_path, baseline_path, force_recalc, comparator_jar = "~/comparator.jar"):
        self.force_recalc = force_recalc
        self.reference = reference_path
        self.comparator_jar = comparator_jar
        self.results = []
        self.__mean_lev_dist = 0
        self.run_calc_stats_task(baseline_path)
        self.baseline = pd.read_csv(Comparator.stats_filename(baseline_path), sep="\t")
        self.baseline.set_index("id", inplace=True)
        self.print_baseline_stats()

    def add(self, path):
        self.run_calc_stats_task(path)
        self.results.append(self.calc_stats(path))

    def save_results(self, path):
        result = pd.DataFrame(self.results)
        result.set_index("path", inplace=True)
        result = result.T
        result.to_latex(path, float_format="%.4f")



class Mapper(object):
    def __init__(self, reference_path, force_remap):
        self.samtools_cmd = "tmap  mapall -f {__reference}".format(__reference=reference_path) \
                            + " -i {__input_type} -o 1 -O 2 -g 3 -s {__output} -o 0 -v stage1  map1 map2 map3 map4"
        self.samtools_view = "samtools view -h {__input_file}"
        self.zcat_view = "zcat {__input_file}"
        self.cat_view = "cat {__input_file}"
        self.force_remap = force_remap

    @staticmethod
    def get_reads_format(path):
        if path.endswith("sam"):
            return "sam"
        elif path.endswith("bam"):
            return "sam"
        elif path.endswith("fasta.gz") or path.endswith("fasta"):
            return "fasta"
        elif path.endswith("fastq"):
            return "fastq"


    @staticmethod
    def is_sam_or_bam(path):
        return path.endswith("sam") or path.endswith("bam")

    @staticmethod
    def is_fasta_gzip(path):
        return path.endswith("fasta.gz")

    @staticmethod
    def is_fasta_or_fastq(path):
        return path.endswith("fasta") or path.endswith("fastq")

    def get_view_cmd(self, path):
        if Mapper.is_sam_or_bam(path):
            return self.samtools_view
        elif Mapper.is_fasta_gzip(path):
            return self.zcat_view
        elif Mapper.is_fasta_or_fastq(path):
            return self.cat_view
        else:
            exit_with_error("Unknown extension for file " + path)

    def map_reads(self, path):
        mapped_path = self.mapped_filename(path)
        if os.path.isfile(mapped_path) and not self.force_remap:
            return
        view_cmd = self.get_view_cmd(path).format(__input_file=path)
        map_cmd = self.samtools_cmd.format(__output=self.mapped_filename(path), __input_type=self.get_reads_format(path))
        cmd = "{__view_cmd} | {__map_cmd}".format(__view_cmd=view_cmd, __map_cmd=map_cmd)
        print("Running tmap command: " + cmd)
        subprocess.call(cmd, shell=True)

    @staticmethod
    def mapped_filename(path):
        return path + ".mapped.sam"

    def run_task(self, path):
        self.map_reads(path)
        return self.mapped_filename(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute reads quality for correction')
    parser.add_argument('--reference', dest='reference', help="Reference fasta file")
    parser.add_argument("--force", help="Force recalc all stats", dest='recalc', action="store_true")
    parser.add_argument("--reads", help="Baseline reads", dest='reads')
    parser.add_argument("--corrected-reads", help="Corrected reads", dest='corrected_reads', nargs="+")
    args = parser.parse_args()

    mapper = Mapper(args.reference, args.recalc)

    mapped_baseline_path = mapper.run_task(args.reads)
    comparator = Comparator(args.reference, mapped_baseline_path, args.recalc)

    for corrected_path in args.corrected_reads:
        mapped_reads_path = mapper.run_task(corrected_path)
        comparator.add(mapped_reads_path)
    comparator.save_results("corrections_quality.tex")


