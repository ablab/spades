__author__ = 'anton'
import os
import sys

sys.path.append("src/spades_pipeline/common")
sys.path.append("src/spades_pipeline/truspades")
sys.path.append("src/spades_pipeline")
import SeqIO
import parallel_launcher
import reference_construction
from barcode_extraction import extract_barcodes, Barcode

import launch_options
import support


def generate_dataset(input_dir):
    sys.stdout.write("Generating truseq dataset from files in directory " + input_dir + "\n")
    if os.path.exists(input_dir) and os.path.isdir(input_dir):
        files = [os.path.abspath(os.path.join(input_dir, file)) for file in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, file))]
        if len(files) == 0:
            sys.stdout.write("Error: Input directory does not contain reads\n")
            sys.exit(1)
    else:
        launch_options.print_usage()
        sys.exit(1)
    return extract_barcodes(files)

def ReadDataset(file):
    sys.stdout.write("Reading dataset from " + file + "\n")
    if os.path.exists(file) and os.path.isfile(file):
        result = []
        f = open(file, "r")
        lines = f.xreadlines()
        for line in lines:
            line = line.strip()
            if line == "":
                continue
            split = line.split(" ")
            id = split[0]
            datasets = []
            for i in range(1, len(split), 2):
                datasets.append([split[i], split[i + 1]])
            result.append(Barcode(datasets, id))
        f.close()
        return result
    else:
        sys.stderr.write("Error: Dataset file does not exist\n" + file + "\n")
        launch_options.print_usage()
        sys.exit(1)

def print_lines(lines, file):
    f = open(file, "w")
    f.write("\n".join([str(line).strip() for line in lines]))
    f.write("\n")
    f.close()

def print_dataset(dataset, file):
    sys.stdout.write("Printing dataset to " + file + "\n")
    print_lines(dataset, file)

def reads_line(libs):
    result = []
    for i in range(len(libs)):
        result.append("--pe{0}-1".format(str(i + 1)))
        result.append(libs[i][0])
        result.append("--pe{0}-2".format(str(i + 1)))
        result.append(libs[i][1])
    return " ".join(result)

#todo: replace with Job class and make continue its parameter

def command_line(barcode, output_dir, params, continue_launch):
#    logfile = os.path.join(output_dir, "logs", barcode.id + ".out")
    if continue_launch and os.path.exists(os.path.join(output_dir, barcode.id,  "params.txt")):
        result = ["./spades.py", "--truseq", "-o", os.path.join(output_dir, barcode.id), "--continue", " ".join(params)]
    else:
       result = ["./spades.py", "--truseq", "-t", "1", "-o", os.path.join(output_dir, barcode.id), reads_line(barcode.libs), " ".join(params)]
#    result = ["./truspades.py", "-o", os.path.join(output_dir, barcode.id), reads_line(barcode.libs), " ".join(params), "\n"]
    return " ".join(result)

def print_commands(commands, options):
    output_file = os.path.join(options.output_dir, "spades_commands.info")
    sys.stdout.write("Printing commands to " + output_file + "\n")
    print_lines(commands, output_file)

def collect_contigs(dataset, output_dir, format):
    output = open(os.path.join(output_dir, "truseq_long_reads." + format), "w")
    for barcode in dataset:
        file = os.path.join(output_dir, barcode.id, "truseq_long_reads." + format)
        if os.path.exists(file):
            contigs = SeqIO.parse(open(file), format)
            for contig in contigs:
                contig.id = barcode.id + "-" + contig.id
                SeqIO.write(contig, output, format)
        else:
            sys.stderr.write("Warning: could not find assembly results for barcode " + str(barcode.id) + "\n")
    output.close()


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

def RunTruSPAdes(dataset, log_dir, options):
    sys.stdout.write("Launching truSPAdes assembly in " + str(options.threads) + " threads\n")
    sys.stdout.write("You can find logs for separate barcodes in " + log_dir + "\n")
    commands = [(barcode.id, command_line(barcode, options.output_dir, options.spades_options, options.continue_launch))
                for barcode in dataset]
    task = parallel_launcher.ExternalCallTask("", "")
    errors = parallel_launcher.run_in_parallel(task, commands, options.threads)
    if errors != 0:
        sys.stderr.write(str(errors) + " barcodes failed to assemble\n")
    collect_contigs(dataset, options.output_dir, "fasta")
    collect_contigs(dataset, options.output_dir, "fastq")
    sys.stdout.write("Assembled virtual long TruSeq reads can be found in " + os.path.join(options.output_dir,
                                                                                           "truseq_long_reads.fasta") + "\n")


def main(argv):
    options = launch_options.Options(argv)
    if options.help:
        launch_options.print_usage()
        sys.exit(0)
    dataset_file = os.path.join(options.output_dir, "dataset.info")
    if options.input_dir is not None:
        dataset = generate_dataset(options.input_dir)
        if dataset is None:
            sys.stderr.write("Error: could not parse dataset from directory " + options.input_dir + "\n")
            sys.exit(1)
        print_dataset(dataset, dataset_file)
        sys.stdout.write("Dataset generated. See result in " + dataset_file)
    else:
        dataset = ReadDataset(options.dataset_file)
        print_dataset(dataset, dataset_file)
    support.ensure_dir_existence(options.output_dir)
    log_dir = os.path.join(options.output_dir, "logs")
    support.ensure_dir_existence(log_dir)
    # if options.print_commands:
    #     verify_exists(options.output_dir)
#         print_commands(commands, options)
    if options.mode == "run_truspades":
        RunTruSPAdes(dataset, log_dir, options)
    elif options.mode == "construct_subreferences":
        reference_construction.ConstructSubreferences(dataset, options.reference, options.output_dir, options.index, options.threads, log = None)
    sys.stdout.write("TruSPAdes launch successfully finished\n")

if __name__ == '__main__':
    main(sys.argv)
