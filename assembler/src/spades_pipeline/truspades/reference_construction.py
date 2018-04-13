#! /usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


# wraps samtools, used to read SAM/BAM
#import pysam3 as pysam
import logging
import os
import shutil
import sys
import itertools
import traceback

import SeqIO
import alignment
import parallel_launcher
import sam_parser
import support

# for fasta/q I/O

class Record:
    def __init__(self, sam_rec):
        self.rname = sam_rec.tname
        self.left = sam_rec.pos
        self.right = sam_rec.pos + sam_rec.alen
        self.cov = sam_rec.alen

    def Coverage(self):
        return 1.0 * self.cov / (self.right - self.left)

    def Len(self):
        return self.right - self.left

    def Join(self, rec):
        self.left = min(self.left, rec.left)
        self.right = max(self.right, rec.right)
        self.cov += rec.cov

    def __str__(self):
        return str(self.rname) + "(" + str(self.Coverage()) + "): [" + str(self.left) + ", " + str(self.right) + "]"
    
    def __cmp__(self, other):
        if other == None:
            return -1
        if self.rname != other.rname:
            if self.rname < other.rname:
                return -1
            else:
                return 1
        if self.left != other.left:
            return self.left - other.left
        return self.right - other.right

def CollectParts(recs, step, mincov, minlen):
    res = []
    cur = None
    for rec in recs:
        if None == cur or cur.rname != rec.rname or cur.right + step < rec.left:
            if cur != None and cur.Len() > minlen and cur.Coverage() > mincov:
                res.append(cur)
            cur = rec
        else:
            cur.Join(rec)
    if None != cur and cur.Len() > minlen and cur.Coverage() > mincov:
        res.append(cur)
    return res

def PrintParts(recs, out):
    for rec in recs:
        out.write(str(rec) + "\n")

def ReadReference(file):
    result = dict()
    for rec in SeqIO.parse_fasta(open(file, "r")):
        result[rec.id] = rec.seq
    return result


def ConstructSubreferenceFromSam(sam_files):
    #todo: make online
    #todo: use config
    recs = []
    for sam_file in sam_files:
        sam = sam_parser.Samfile(sam_file)
        for rec in sam:
            if rec.pos != -1:
                recs.append(Record(rec))
    recs.sort()
    filtered_recs = CollectParts(recs, 0, 2, 500)
    subreferences = CollectParts(filtered_recs, 20000, 0, 7000)
    return filtered_recs, subreferences


def PrintResults(recs, reference, references_file, coordinates_file):
    aln = open(coordinates_file, "w")
    fasta = open(references_file, "w")
    for rec in recs:
        aln.write(str(rec) + "\n")
        sequence = reference[rec.rname][rec.left:rec.right]
        rec_id = str(rec.rname) + "_(" + str(rec.left) + "-" + str(rec.right)+")"
        SeqIO.write(SeqIO.SeqRecord(sequence, rec_id), fasta, "fasta")
    aln.close()
    fasta.close()

def PrintAll(subreferences, reference, output_dir):
    references_dir = os.path.join(output_dir, "references")
    coordinates_dir = os.path.join(output_dir, "coordinates")
    os.mkdir(references_dir)
    os.mkdir(coordinates_dir)
    for id, subreference in subreferences:
        if subreference != None:
            PrintResults(subreference, reference, os.path.join(references_dir, str(id) + ".fasta"), os.path.join(coordinates_dir, str(id) + ".aln"))

def AlignToReference(datasets, sam_dir, bwa_command, log, index, threads = 1):
    if os.path.exists(sam_dir):
        shutil.rmtree(sam_dir)
    os.makedirs(sam_dir)
    for dataset_id, reads in datasets:
        log.info("\nProcessing barcode " + str(dataset_id))
        dataset_work_dir = os.path.join(sam_dir, str(dataset_id))
        sam_files = alignment.align_bwa_pe_libs(bwa_command, index, reads, dataset_work_dir, log, threads)
        yield (dataset_id, sam_files)

def ReadDataset(dataset_file):
    dataset_file = dataset_file.split(":")
    dataset_lines = [line.strip().split(" ") for line in open(dataset_file[0], "r").readlines() if line.strip()]
    datasets = [(line[0], zip(line[1::2], line[2::2])) for line in dataset_lines]
    if len(dataset_file) > 1:
        datasets = filter(lambda x, y: x.startswith(dataset_file[1]), datasets)
    return datasets

def ConstructSubreferences(datasets, reference_file, output_dir, index = None, threads = 1, log = None):
    bwa_command = "bin/spades-bwa"
    if log == None:
        log = logging.getLogger('reference_construction')
        log.setLevel(logging.INFO)
        console = logging.StreamHandler(sys.stderr)
        console.setFormatter(logging.Formatter('%(message)s'))
        console.setLevel(logging.INFO)
        log.addHandler(console)
    support.ensure_dir_existence(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if index == None:
        log.info("Constructing index\n")
        index = alignment.index_bwa(bwa_command, log, reference_file, os.path.join(output_dir, "bwa_index"), "bwtsw")
    sam_dir = os.path.join(output_dir, "alignments")
    log.info("Aligning barcodes\n")
    sam_files = AlignToReference(datasets, sam_dir, bwa_command, log, index, threads)
    subreference_dir = os.path.join(output_dir, "subreferences")
    filtered_dir = os.path.join(output_dir, "filtered")
    support.recreate_dir(subreference_dir)
    support.recreate_dir(filtered_dir)
    log.info("Constructing subreferences")
    subreferences_list = [(barcode_id, ConstructSubreferenceFromSam(barcode_sam)) for barcode_id, barcode_sam in sam_files]
    log.info("Reading reference")
    reference = ReadReference(reference_file)
    log.info("Printing output")
    PrintAll([(barcode, filtered) for barcode, (filtered, subreference) in subreferences_list], reference, filtered_dir)
    PrintAll([(barcode, subreference) for barcode, (filtered, subreference) in subreferences_list], reference, subreference_dir)
    log.info("Subreference construction finished. See results in " + output_dir)

if __name__ == '__main__':
    # ConstructSubreference(sys.argv[1], "r", ReadReference(sys.argv[2]), sys.argv[3])
    ConstructSubreferences(ReadDataset(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]))
