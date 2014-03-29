#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import support

SUPPORTED_PYTHON_VERSIONS = ['2.4', '2.5', '2.6', '2.7', '3.2', '3.3']
# allowed reads extensions for BayesHammer and for thw whole SPAdes pipeline
BH_ALLOWED_READS_EXTENSIONS = ['.fq', '.fastq']
CONTIGS_ALLOWED_READS_EXTENSIONS = ['.fa', '.fasta']
ALLOWED_READS_EXTENSIONS = BH_ALLOWED_READS_EXTENSIONS + CONTIGS_ALLOWED_READS_EXTENSIONS
# reads could be gzipped
BH_ALLOWED_READS_EXTENSIONS += [x + '.gz' for x in BH_ALLOWED_READS_EXTENSIONS]
CONTIGS_ALLOWED_READS_EXTENSIONS += [x + '.gz' for x in CONTIGS_ALLOWED_READS_EXTENSIONS]
ALLOWED_READS_EXTENSIONS += [x + '.gz' for x in ALLOWED_READS_EXTENSIONS]

# we support up to MAX_LIBS_NUMBER paired-end libs and MAX_LIBS_NUMBER mate-pair libs
MAX_LIBS_NUMBER = 5
# other libs types:
LONG_READS_TYPES = ["pacbio", "sanger", "trusted-contigs", "untrusted-contigs"]

#other constants
MIN_K = 1
MAX_K = 127
THRESHOLD_FOR_BREAKING_SCAFFOLDS = 3
THRESHOLD_FOR_BREAKING_ADDITIONAL_CONTIGS = 10

#default values constants
THREADS = 16
MEMORY = 250
K_MERS_SHORT = [21,33,55]
K_MERS_150 = [21,33,55,77]
K_MERS_250 = [21,33,55,77,99,127]
ITERATIONS = 1
TMP_DIR = "tmp"

### START OF OPTIONS
# basic options
output_dir = None
single_cell = False
iontorrent = False

# pipeline options
only_error_correction = False
only_assembler = False
disable_gzip_output = None
disable_rr = None
careful = None
diploid_mode = False

# advanced options
continue_mode = False
developer_mode = None
dataset_yaml_filename = None
threads = None
memory = None
tmp_dir = None
k_mers = None
qvoffset = None # auto-detect by default

# hidden options
mismatch_corrector = None
reference = None
iterations = None
bh_heap_check = None
spades_heap_check = None
read_buffer_size = None
### END OF OPTIONS

# for restarting SPAdes
restart_from = None
restart_careful = None
restart_mismatch_corrector = None
restart_disable_gzip_output = None
restart_disable_rr = None
restart_threads = None
restart_memory = None
restart_tmp_dir = None
restart_k_mers = None
restart_qvoffset = None
restart_developer_mode = None
restart_reference = None
restart_read_buffer_size = None

dict_of_prefixes = dict()

# list of spades.py options
long_options = "12= threads= memory= tmp-dir= iterations= phred-offset= sc iontorrent "\
               "only-error-correction only-assembler "\
               "disable-gzip-output disable-gzip-output:false disable-rr disable-rr:false " \
               "help test debug debug:false reference= config-file= dataset= "\
               "bh-heap-check= spades-heap-check= read-buffer-size= help-hidden "\
               "mismatch-correction mismatch-correction:false careful careful:false "\
               "continue restart-from= diploid".split()
short_options = "o:1:2:s:k:t:m:i:h"

# adding multiple paired-end, mate-pair and other (long reads) libraries support
reads_options = []
for i in range(MAX_LIBS_NUMBER):
    for type in ["pe", "mp"]:
        reads_options += ("%s%d-1= %s%d-2= %s%d-12= %s%d-s= %s%d-rf %s%d-fr %s%d-ff" % tuple([type, i + 1] * 7)).split()
reads_options += list(map(lambda x: x + '=', LONG_READS_TYPES))
long_options += reads_options
# for checking whether option corresponds to reads or not
reads_options = list(map(lambda x:"--" + x.split('=')[0], reads_options))
reads_options += ["--12", "-1", "-2", "-s"]


def usage(spades_version, show_hidden=False, dipspades=False):
    if not dipspades:
        sys.stderr.write("SPAdes genome assembler v." + str(spades_version) + "\n")
    else:
        sys.stderr.write("dipSPAdes 1.0: genome assembler designed for diploid genomes with high heterozygosity rate\n\n")
    sys.stderr.write("Usage: " + str(sys.argv[0]) + " [options] -o <output_dir>" + "\n")
    sys.stderr.write("" + "\n")
    sys.stderr.write("Basic options:" + "\n")
    sys.stderr.write("-o\t<output_dir>\tdirectory to store all the resulting files (required)" + "\n")
    if not dipspades:
        sys.stderr.write("--sc\t\t\tthis flag is required for MDA (single-cell) data" + "\n")
    sys.stderr.write("--iontorrent\t\tthis flag is required for IonTorrent data" + "\n")
    sys.stderr.write("--test\t\t\truns SPAdes on toy dataset" + "\n")
    sys.stderr.write("-h/--help\t\tprints this usage message" + "\n")

    sys.stderr.write("" + "\n")
    if not dipspades:
        sys.stderr.write("Input data:" + "\n")
    else:
        sys.stderr.write("Input reads:" + "\n")
    sys.stderr.write("--12\t<filename>\tfile with interlaced forward and reverse"\
                         " paired-end reads" + "\n")
    sys.stderr.write("-1\t<filename>\tfile with forward paired-end reads" + "\n")
    sys.stderr.write("-2\t<filename>\tfile with reverse paired-end reads" + "\n")
    sys.stderr.write("-s\t<filename>\tfile with unpaired reads" + "\n")
    sys.stderr.write("--pe<#>-12\t<filename>\tfile with interlaced"\
                         " reads for paired-end library number <#> (<#> = 1,2,3,4,5)" + "\n")
    sys.stderr.write("--pe<#>-1\t<filename>\tfile with forward reads"\
                         " for paired-end library number <#> (<#> = 1,2,3,4,5)" + "\n")
    sys.stderr.write("--pe<#>-2\t<filename>\tfile with reverse reads"\
                         " for paired-end library number <#> (<#> = 1,2,3,4,5)" + "\n")
    sys.stderr.write("--pe<#>-s\t<filename>\tfile with unpaired reads"\
                         " for paired-end library number <#> (<#> = 1,2,3,4,5)" + "\n")
    sys.stderr.write("--pe<#>-<or>\torientation of reads"\
                         " for paired-end library number <#> (<#> = 1,2,3,4,5; <or> = fr, rf, ff)" + "\n")
    sys.stderr.write("--mp<#>-12\t<filename>\tfile with interlaced"\
                         " reads for mate-pair library number <#> (<#> = 1,2,3,4,5)" + "\n")
    sys.stderr.write("--mp<#>-1\t<filename>\tfile with forward reads"\
                         " for mate-pair library number <#> (<#> = 1,2,3,4,5)" + "\n")
    sys.stderr.write("--mp<#>-2\t<filename>\tfile with reverse reads"\
                         " for mate-pair library number <#> (<#> = 1,2,3,4,5)" + "\n")
    sys.stderr.write("--mp<#>-s\t<filename>\tfile with unpaired reads"\
                         " for mate-pair library number <#> (<#> = 1,2,3,4,5)" + "\n")
    sys.stderr.write("--mp<#>-<or>\torientation of reads"\
                         " for mate-pair library number <#> (<#> = 1,2,3,4,5; <or> = fr, rf, ff)" + "\n")
    sys.stderr.write("--sanger\t<filename>\tfile with Sanger reads\n")
    sys.stderr.write("--pacbio\t<filename>\tfile with PacBio reads\n")
    sys.stderr.write("--trusted-contigs\t<filename>\tfile with trusted contigs\n")
    sys.stderr.write("--untrusted-contigs\t<filename>\tfile with untrusted contigs\n")
    if dipspades:
        sys.stderr.write("Input haplocontigs:" + "\n")
        sys.stderr.write("--hap\t<filename>\tfile with haplocontigs" + "\n")

    sys.stderr.write("" + "\n")
    sys.stderr.write("Pipeline options:" + "\n")
    if not dipspades:
        sys.stderr.write("--only-error-correction\truns only read error correction"\
                             " (without assembling)" + "\n")
    sys.stderr.write("--only-assembler\truns only assembling (without read error"\
                         " correction)" + "\n")
    if not dipspades:
        sys.stderr.write("--careful\t\ttries to reduce number"\
                             " of mismatches and short indels" + "\n")
        sys.stderr.write("--continue\t\tcontinue run from the last available check-point" + "\n")
        sys.stderr.write("--restart-from\t<cp>\trestart run with updated options and from the specified check-point ('ec', 'as', 'k<int>', 'mc')" + "\n")
    sys.stderr.write("--disable-gzip-output\tforces error correction not to"\
                         " compress the corrected reads" + "\n")
    sys.stderr.write("--disable-rr\t\tdisables repeat resolution stage"\
                     " of assembling" + "\n")

    if dipspades:
        sys.stderr.write("" + "\n")
        sys.stderr.write("DipSPAdes options:" + "\n")
        sys.stderr.write("--expect-gaps\tindicate that significant number of gaps in coverage is expected" + "\n")
        sys.stderr.write("--expect-rearrangements\tindicate that significant number of rearrangngements between haplomes of diploid genome is expected" + "\n")

    sys.stderr.write("" + "\n")
    sys.stderr.write("Advanced options:" + "\n")
    sys.stderr.write("--dataset\t<filename>\tfile with dataset description in YAML format" + "\n")
    sys.stderr.write("-t/--threads\t<int>\t\tnumber of threads" + "\n")
    sys.stderr.write("\t\t\t\t[default: %s]\n" % THREADS)
    sys.stderr.write("-m/--memory\t<int>\t\tRAM limit for SPAdes in Gb"\
                         " (terminates if exceeded)" + "\n")
    sys.stderr.write("\t\t\t\t[default: %s]\n" % MEMORY)
    sys.stderr.write("--tmp-dir\t<dirname>\tdirectory for temporary files" + "\n")
    sys.stderr.write("\t\t\t\t[default: <output_dir>/tmp]" + "\n")
    sys.stderr.write("-k\t\t<int,int,...>\tcomma-separated list of k-mer sizes"\
                         " (must be odd and" + "\n")
    sys.stderr.write("\t\t\t\tless than " + str(MAX_K + 1) + ") [default: 'auto']" + "\n") # ",".join(map(str, k_mers_short))
    sys.stderr.write("--phred-offset\t<33 or 64>\tPHRED quality offset in the"\
                         " input reads (33 or 64)" + "\n")
    sys.stderr.write("\t\t\t\t[default: auto-detect]" + "\n")    

    if show_hidden:
        sys.stderr.write("" + "\n")
        sys.stderr.write("HIDDEN options:" + "\n")
        sys.stderr.write("--debug\t\t\t\truns SPAdes in debug mode (keeps intermediate output)" + "\n")
        sys.stderr.write("--mismatch-correction\t\truns post processing correction"\
                             " of mismatches and short indels" + "\n")
        sys.stderr.write("--reference\t<filename>\tfile with reference for deep analysis"\
                             " (only in debug mode)" + "\n")
        sys.stderr.write("-i/--iterations\t<int>\t\tnumber of iterations for read error"\
                             " correction [default: %s]\n" % ITERATIONS)
        sys.stderr.write("--read-buffer-size\t<int>\t\tsets size of read buffer for graph construction")
        sys.stderr.write("--bh-heap-check\t\t<value>\tsets HEAPCHECK environment variable"\
                             " for BayesHammer" + "\n")
        sys.stderr.write("--spades-heap-check\t<value>\tsets HEAPCHECK environment variable"\
                             " for SPAdes" + "\n")
        sys.stderr.write("--help-hidden\tprints this usage message with all hidden options" + "\n")

    sys.stderr.flush()


def auto_K_allowed():
    return not k_mers and not single_cell and not iontorrent # kmers were set by default, not SC, and not IonTorrent data


def set_default_values():
    global threads
    global memory
    global iterations
    global disable_gzip_output
    global disable_rr
    global careful
    global mismatch_corrector
    global developer_mode
    global qvoffset
    global tmp_dir

    if threads is None:
        threads = THREADS
    if memory is None:
        memory = MEMORY
    if iterations is None:
        iterations = ITERATIONS
    if disable_gzip_output is None:
        disable_gzip_output = False
    if disable_rr is None:
        disable_rr = False
    if careful is None:
        careful = False
    if mismatch_corrector is None:
        mismatch_corrector = False
    if developer_mode is None:
        developer_mode = False
    if qvoffset == 'auto':
        qvoffset = None
    if tmp_dir is None:
        tmp_dir = os.path.join(output_dir, TMP_DIR)


def set_test_options():
    global output_dir
    global single_cell

    output_dir = os.path.abspath('spades_test')
    single_cell = False


def save_restart_options(log):
    if dataset_yaml_filename:
        support.error("you cannot specify --dataset with --restart-from option!", log)
    if single_cell:
        support.error("you cannot specify --sc with --restart-from option!", log)
    if iontorrent:
        support.error("you cannot specify --iontorrent with --restart-from option!", log)
    if only_assembler:
        support.error("you cannot specify --only-assembler with --restart-from option!", log)
    if only_error_correction:
        support.error("you cannot specify --only-error-correction with --restart-from option!", log)

    global restart_k_mers
    global restart_careful
    global restart_mismatch_corrector
    global restart_disable_gzip_output
    global restart_disable_rr
    global restart_threads
    global restart_memory
    global restart_tmp_dir
    global restart_qvoffset
    global restart_developer_mode
    global restart_reference
    global restart_read_buffer_size

    restart_k_mers = k_mers
    restart_careful = careful
    restart_mismatch_corrector = mismatch_corrector
    restart_disable_gzip_output = disable_gzip_output
    restart_disable_rr = disable_rr
    restart_threads = threads
    restart_memory = memory
    restart_tmp_dir = tmp_dir
    restart_qvoffset = qvoffset
    restart_developer_mode = developer_mode
    restart_reference = reference
    restart_read_buffer_size = read_buffer_size


def load_restart_options():
    global k_mers
    global careful
    global mismatch_corrector
    global disable_gzip_output
    global disable_rr
    global threads
    global memory
    global tmp_dir
    global qvoffset
    global developer_mode
    global reference
    global read_buffer_size

    if restart_k_mers:
        if restart_k_mers == 'auto':
            k_mers = None  # set by default
        else:
            k_mers = restart_k_mers
    if restart_careful is not None:
        careful = restart_careful
    if restart_mismatch_corrector is not None:
        mismatch_corrector = restart_mismatch_corrector
    if disable_gzip_output is not None:
        disable_gzip_output = restart_disable_gzip_output
    if restart_disable_rr is not None:
        disable_rr = restart_disable_rr
    if restart_threads is not None:
        threads = restart_threads
    if restart_memory is not None:
        memory = restart_memory
    if restart_tmp_dir is not None:
        tmp_dir = restart_tmp_dir
    if restart_qvoffset is not None:
        qvoffset = restart_qvoffset
    if restart_developer_mode is not None:
        developer_mode = restart_developer_mode
    if restart_reference is not None:
        reference = restart_reference
    if restart_read_buffer_size is not None:
        read_buffer_size = restart_read_buffer_size
