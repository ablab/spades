#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################
import gzip
import logging

import os
import sys
import shutil
import spades_init
spades_init.init()
truspades_home = spades_init.spades_home
spades_home = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
spades_version = spades_init.spades_version

import SeqIO  # TODO: add to ext/scr/python_libs
import parallel_launcher
import reference_construction
import launch_options
import support
import barcode_extraction

def generate_dataset(input_dirs, log):
    log.info("Generating truseq dataset from input directories:")
    for idir in input_dirs:
        log.info("\t" + idir)
    for input_dir in input_dirs:
        if not os.path.exists(input_dir) or not os.path.isdir(input_dir):
            log.info("Input directory " + input_dir + " does not exist")
            sys.exit(1)
        files = [os.path.abspath(os.path.join(input_dir, file)) for file in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, file))]
        if len(files) == 0:
            log.info("Error: Input directory does not contain reads")
            sys.exit(1)
    return barcode_extraction.ExtractBarcodes(input_dirs)

def reads_line(libs):
    result = []
    for i in range(len(libs)):
        result.append("-1")
        result.append(libs[i][0])
        result.append("-2")
        result.append(libs[i][1])
    return " ".join(result)

#todo: replace with Job class and make continue its parameter

def command_line(barcode, output_dir, params, continue_launch):
#    logfile = os.path.join(output_dir, "logs", barcode.id + ".out")
    if continue_launch and os.path.exists(os.path.join(output_dir, barcode.id,  "params.txt")):
        result = ["python " + os.path.join(spades_home, "spades.py"), "--truseq", "-o", os.path.join(output_dir, barcode.id), "--continue", params]
    else:
       result = ["python " + os.path.join(spades_home, "spades.py"), "--truseq", "-t", "1", "-o", os.path.join(output_dir, barcode.id), reads_line(barcode.libs), params]
#    result = ["./truspades.py", "-o", os.path.join(output_dir, barcode.id), reads_line(barcode.libs), " ".join(params), "\n"]
    return " ".join(result)

def print_commands(commands, options, log):
    output_file = os.path.join(options.output_dir, "spades_commands.info")
    log.info("Printing commands to " + output_file)
    open(output_file, "w").write("\n".join([str(line).strip() for line in commands]) + "\n")

def collect_contigs(dataset, barcodes_dir, output_base, format):
    output = open(output_base + "." + format, "w")
    for barcode in dataset:
        file = os.path.join(barcodes_dir, barcode.id, "truseq_long_reads." + format)
        if os.path.exists(file):
            contigs = SeqIO.parse(open(file), format)
            for contig in contigs:
                contig.id = barcode.id + "-" + contig.id
                SeqIO.write(contig, output, format)
    output.close()

def check_results(dataset, output_dir, log):
    for barcode in dataset:
        if not os.path.exists(os.path.join(output_dir, barcode.id, "truseq_long_reads.fastq")):
            log.info("Warning: could not find assembly results for barcode " + str(barcode.id))

def bwa_command_line(barcode, output_dir, index, threads):
    return " ".join(["bwa", "mem", "-t", str(threads), index, barcode.libs[0][0], barcode.libs[0][1]])

class ReferenceConstructionLauncher:
    def __init__(self, reference,  sam_dir, result_dir):
        self.reference = reference
        self.sam_dir = sam_dir
        self.result_dir = result_dir
    def __call__(self, barcode_id):
       reference_construction.CounstructSubreference(os.path.join(self.sam_dir, barcode_id + ".sam"), reference_construction.ReadReference(self.reference), os.path.join(self.result_dir, barcode_id))
       return 0
        
def ConstructSubreferences(dataset, options):
    reference_construction.ConstructSubreferences(dataset, options.reference, options.output_dir, options.index, options.threads, log = None)

def RunTruSPAdes(dataset, log_dir, options, log):
    log.info("Launching truSPAdes assembly in " + str(options.threads) + " threads")
    log.info("You can find logs for separate barcodes in " + log_dir)
    barcodes_dir = os.path.join(options.output_dir, "barcodes")
    support.ensure_dir_existence(barcodes_dir)
    commands = [(barcode.id, command_line(barcode, barcodes_dir, options.spades_options, options.continue_launch))
                for barcode in dataset]
    task = parallel_launcher.ExternalCallTask(os.path.join(log_dir, "{0}.log"), "", log.name)
    errors = parallel_launcher.run_in_parallel(task, commands, options.threads)
    if errors != 0:
        log.info(str(errors) + " barcodes failed to assemble")
    check_results(dataset, barcodes_dir, log)
    output_base = os.path.join(options.output_dir, "TSLR")
    collect_contigs(dataset, barcodes_dir, output_base, "fasta")
    collect_contigs(dataset, barcodes_dir, output_base, "fastq")
    log.info("Assembled virtual long TruSeq reads can be found in " + os.path.join(options.output_dir,
                                                                                           "TSLR.fasta"))
    if options.clean:
        SaveContigs(barcodes_dir, dataset, "fasta")
        SaveContigs(barcodes_dir, dataset, "fastq")
        for barcode in dataset:
            shutil.rmtree(os.path.join(barcodes_dir, barcode.id))


def SaveContigs(barcodes_dir, dataset, format):
    contig_dir = os.path.join(barcodes_dir, format)
    support.ensure_dir_existence(contig_dir)
    for barcode in dataset:
        if os.path.isfile(os.path.join(barcodes_dir, barcode.id, "truseq_long_reads." + format)):
            shutil.copyfileobj(open(os.path.join(barcodes_dir, barcode.id, "truseq_long_reads." + format), "rb"), gzip.open(os.path.join(contig_dir, barcode.id + "." + format + ".gz"), "wb"))


def create_log(options):
    log = logging.getLogger('truspades')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)
    log_filename = os.path.join(options.output_dir, "truspades.log")
    if options.continue_launch:
        log_handler = logging.FileHandler(log_filename, mode='a')
    else:
        log_handler = logging.FileHandler(log_filename, mode='w')
    log.addHandler(log_handler)
    return log


def CheckTestSuccess(options, log):
    output = os.path.join(options.output_dir, "TSLR.fasta")
    if not os.path.isfile(output):
        log.info("TruSPAdes test launch failed: can not find output files")
        sys.exit(1)
    if not (os.path.getsize(output) > 20000 and os.path.getsize(output) < 20100):
        log.info("TruSPAdes test launch failed: incorrect output files")
        sys.exit(1)
    log.info("TruSPAdes test passed correctly")



def main(argv):
    options = launch_options.Options(argv, spades_home, truspades_home, spades_version)
    support.ensure_dir_existence(options.output_dir)
    if options.test and not options.continue_launch:
        support.recreate_dir(options.output_dir)
    log = create_log(options)
    dataset_file = os.path.join(options.output_dir, "dataset.info")
    if options.continue_launch:
        dataset = barcode_extraction.ReadDataset(dataset_file, log)
    elif options.input_dirs is not None:
        dataset = generate_dataset(options.input_dirs, log)
        if dataset is None:
            log.info("Error: could not parse dataset from input directories\n")
            sys.exit(1)
        barcode_extraction.print_dataset(dataset, dataset_file, log)
        log.info("Dataset generated. See result in " + dataset_file)
    else:
        dataset = barcode_extraction.ReadDataset(options.dataset_file, log)
        barcode_extraction.print_dataset(dataset, dataset_file, log)
    log_dir = os.path.join(options.output_dir, "logs")
    support.ensure_dir_existence(log_dir)
    # if options.print_commands:
    #     verify_exists(options.output_dir)
#         print_commands(commands, options)
    if options.mode == "run_truspades":
        RunTruSPAdes(dataset, log_dir, options, log)
    elif options.mode == "construct_subreferences":
        reference_construction.ConstructSubreferences(dataset, options.reference, options.output_dir, options.index, options.threads, log = None)
    log.info("TruSPAdes launch successfully finished")
    if options.test:
        CheckTestSuccess(options, log)

if __name__ == '__main__':
    main(sys.argv)
