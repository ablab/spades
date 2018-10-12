#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import support
from os.path import basename

SUPPORTED_PYTHON_VERSIONS = ['2.4-2.7', '3.2+']  # major.minor format only, close ('-') and open ('+') ranges allowed
# allowed reads extensions for BayesHammer and for thw whole SPAdes pipeline
BH_ALLOWED_READS_EXTENSIONS = ['.fq', '.fastq', '.bam']
CONTIGS_ALLOWED_READS_EXTENSIONS = ['.fa', '.fasta']
ALLOWED_READS_EXTENSIONS = BH_ALLOWED_READS_EXTENSIONS + CONTIGS_ALLOWED_READS_EXTENSIONS
# reads could be gzipped
BH_ALLOWED_READS_EXTENSIONS += [x + '.gz' for x in BH_ALLOWED_READS_EXTENSIONS]
CONTIGS_ALLOWED_READS_EXTENSIONS += [x + '.gz' for x in CONTIGS_ALLOWED_READS_EXTENSIONS]
ALLOWED_READS_EXTENSIONS += [x + '.gz' for x in ALLOWED_READS_EXTENSIONS]

# we support up to MAX_LIBS_NUMBER libs for each type of short-reads libs
MAX_LIBS_NUMBER = 9
OLD_STYLE_READS_OPTIONS = ["--12", "-1", "-2", "-s", "--merged"]
SHORT_READS_TYPES = {"pe": "paired-end", "s": "single", "mp": "mate-pairs", "hqmp": "hq-mate-pairs", "nxmate": "nxmate"}
# other libs types:
LONG_READS_TYPES = ["pacbio", "sanger", "nanopore", "tslr", "trusted-contigs", "untrusted-contigs"]

# final contigs and scaffolds names
contigs_name = "contigs.fasta"
scaffolds_name = "scaffolds.fasta"
assembly_graph_name = "assembly_graph.fastg"
assembly_graph_name_gfa = "assembly_graph_with_scaffolds.gfa"
contigs_paths = "contigs.paths"
scaffolds_paths = "scaffolds.paths"
transcripts_name = "transcripts.fasta"
transcripts_paths = "transcripts.paths"
filtering_types = ["hard", "soft", "default"]

#other constants
MIN_K = 1
MAX_K = 127
RNA_MIN_K = 29
RNA_MAX_LOWER_K = 55
THRESHOLD_FOR_BREAKING_SCAFFOLDS = 3
THRESHOLD_FOR_BREAKING_ADDITIONAL_CONTIGS = 10

#default values constants
THREADS = 16
MEMORY = 250
K_MERS_RNA = [33,49]
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
meta = False
rna = False
large_genome = False
test_mode = False
plasmid = False

# pipeline options
only_error_correction = False
only_assembler = False
disable_gzip_output = None
disable_rr = None
careful = None

# advanced options
continue_mode = False
developer_mode = None
dataset_yaml_filename = None
threads = None
memory = None
tmp_dir = None
k_mers = None
qvoffset = None  # auto-detect by default
cov_cutoff = 'off'  # default is 'off'

# hidden options
save_gp = False
mismatch_corrector = None
reference = None
series_analysis = None
configs_dir = None
iterations = None
bh_heap_check = None
spades_heap_check = None
read_buffer_size = None
lcer_cutoff = None
read_cov_threshold = None
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
original_k_mers = None
restart_qvoffset = None
restart_cov_cutoff = None
restart_developer_mode = None
restart_reference = None
restart_configs_dir = None
restart_read_buffer_size = None
restart_fast = None

# for running to specific check-point
stop_after = None
run_completed = False

#truseq options
truseq_mode = False
correct_scaffolds = False
run_truseq_postprocessing = False

#rna options
strand_specificity = None  # None, 'rf', 'fr' are possible
fast = None

dict_of_prefixes = dict()
dict_of_rel2abs = dict()

# list of spades.py options
long_options = "12= merged= threads= memory= tmp-dir= iterations= phred-offset= sc iontorrent meta large-genome rna plasmid "\
               "ss-fr ss-rf fast fast:false "\
               "only-error-correction only-assembler "\
               "disable-gzip-output disable-gzip-output:false disable-rr disable-rr:false " \
               "help version test debug debug:false reference= series-analysis= config-file= dataset= "\
               "bh-heap-check= spades-heap-check= read-buffer-size= help-hidden "\
               "mismatch-correction mismatch-correction:false careful careful:false save-gp save-gp:false "\
               "continue restart-from= truseq cov-cutoff= hidden-cov-cutoff= read-cov-threshold= " \
               "configs-dir= stop-after=".split()
short_options = "o:1:2:s:k:t:m:i:hv"

# adding multiple paired-end, mate-pair and other (long reads) libraries support
reads_options = []
for i in range(MAX_LIBS_NUMBER):
    for type in SHORT_READS_TYPES.keys():
        if type == 's':  # single
            reads_options += ["s%d=" % (i+1)]
        elif type == 'nxmate':  # special case: only left and right reads
            reads_options += ("%s%d-1= %s%d-2=" % tuple([type, i + 1] * 2)).split()
        else:  # paired-end, mate-pairs, hq-mate-pairs
            reads_options += ("%s%d-1= %s%d-2= %s%d-12= %s%d-s= %s%d-rf %s%d-fr %s%d-ff" % tuple([type, i + 1] * 7)).split()
            if type == 'pe':  # special case: paired-end may include merged reads (-m)
                reads_options += ["%s%d-m=" % (type, i+1)]
reads_options += list(map(lambda x: x + '=', LONG_READS_TYPES))
long_options += reads_options
# for checking whether option corresponds to reads or not
reads_options = list(map(lambda x: "--" + x.split('=')[0], reads_options))
reads_options += OLD_STYLE_READS_OPTIONS


def get_mode():
    mode = None
    if basename(sys.argv[0]) == "rnaspades.py":
        mode = 'rna'
    elif basename(sys.argv[0]) == "plasmidspades.py":
        mode = 'plasmid'
    elif basename(sys.argv[0]) == "metaspades.py":
        mode = 'meta'
    return mode


def version(spades_version, mode=None):
    sys.stdout.write("SPAdes v" + str(spades_version))
    if mode is None:
        mode = get_mode()
    if mode is not None:
        sys.stdout.write(" [" + mode + "SPAdes mode]")
    sys.stdout.write("\n")
    sys.stdout.flush()


def usage(spades_version, show_hidden=False, mode=None):
    sys.stderr.write("SPAdes genome assembler v" + str(spades_version))
    if mode is None:
        mode = get_mode()
    if mode is not None:
        sys.stderr.write(" [" + mode + "SPAdes mode]")
    sys.stderr.write("\n\n")
    sys.stderr.write("Usage: " + str(sys.argv[0]) + " [options] -o <output_dir>" + "\n")
    sys.stderr.write("" + "\n")
    sys.stderr.write("Basic options:" + "\n")
    sys.stderr.write("-o\t<output_dir>\tdirectory to store all the resulting files (required)" + "\n")
    if mode is None:  # nothing special, just regular spades.py
        sys.stderr.write("--sc\t\t\tthis flag is required for MDA (single-cell) data" + "\n")
        sys.stderr.write("--meta\t\t\tthis flag is required for metagenomic sample data" + "\n")
        sys.stderr.write("--rna\t\t\tthis flag is required for RNA-Seq data \n")
        sys.stderr.write("--plasmid\t\truns plasmidSPAdes pipeline for plasmid detection \n")

    sys.stderr.write("--iontorrent\t\tthis flag is required for IonTorrent data" + "\n")
    sys.stderr.write("--test\t\t\truns SPAdes on toy dataset" + "\n")
    sys.stderr.write("-h/--help\t\tprints this usage message" + "\n")
    sys.stderr.write("-v/--version\t\tprints version" + "\n")

    sys.stderr.write("" + "\n")
    sys.stderr.write("Input data:" + "\n")

    sys.stderr.write("--12\t<filename>\tfile with interlaced forward and reverse"\
                         " paired-end reads" + "\n")
    sys.stderr.write("-1\t<filename>\tfile with forward paired-end reads" + "\n")
    sys.stderr.write("-2\t<filename>\tfile with reverse paired-end reads" + "\n")
    sys.stderr.write("-s\t<filename>\tfile with unpaired reads" + "\n")
    sys.stderr.write("--merged\t<filename>\tfile with merged forward and reverse paired-end reads" + "\n")
    if mode == "meta":
        allowed_lib_ids = "1"
    else:
        allowed_lib_ids = "1,2,...," + str(MAX_LIBS_NUMBER)
    sys.stderr.write("--pe<#>-12\t<filename>\tfile with interlaced"\
                         " reads for paired-end library number <#> (<#> = %s)" % allowed_lib_ids + "\n")
    sys.stderr.write("--pe<#>-1\t<filename>\tfile with forward reads"\
                         " for paired-end library number <#> (<#> = %s)" % allowed_lib_ids + "\n")
    sys.stderr.write("--pe<#>-2\t<filename>\tfile with reverse reads"\
                         " for paired-end library number <#> (<#> = %s)" % allowed_lib_ids + "\n")
    sys.stderr.write("--pe<#>-s\t<filename>\tfile with unpaired reads"\
                         " for paired-end library number <#> (<#> = %s)" % allowed_lib_ids + "\n")
    sys.stderr.write("--pe<#>-m\t<filename>\tfile with merged reads"\
                         " for paired-end library number <#> (<#> = %s)" % allowed_lib_ids + "\n")
    sys.stderr.write("--pe<#>-<or>\torientation of reads"\
                         " for paired-end library number <#> (<#> = %s; <or> = fr, rf, ff)" % allowed_lib_ids + "\n")
    sys.stderr.write("--s<#>\t\t<filename>\tfile with unpaired reads"\
                     " for single reads library number <#> (<#> = %s)" % allowed_lib_ids + "\n")
    if mode not in ["rna", "meta"]:
        sys.stderr.write("--mp<#>-12\t<filename>\tfile with interlaced"\
                             " reads for mate-pair library number <#> (<#> = 1,2,..,9)" + "\n")
        sys.stderr.write("--mp<#>-1\t<filename>\tfile with forward reads"\
                             " for mate-pair library number <#> (<#> = 1,2,..,9)" + "\n")
        sys.stderr.write("--mp<#>-2\t<filename>\tfile with reverse reads"\
                             " for mate-pair library number <#> (<#> = 1,2,..,9)" + "\n")
        sys.stderr.write("--mp<#>-s\t<filename>\tfile with unpaired reads"\
                             " for mate-pair library number <#> (<#> = 1,2,..,9)" + "\n")
        sys.stderr.write("--mp<#>-<or>\torientation of reads"\
                             " for mate-pair library number <#> (<#> = 1,2,..,9; <or> = fr, rf, ff)" + "\n")
        sys.stderr.write("--hqmp<#>-12\t<filename>\tfile with interlaced"\
                         " reads for high-quality mate-pair library number <#> (<#> = 1,2,..,9)" + "\n")
        sys.stderr.write("--hqmp<#>-1\t<filename>\tfile with forward reads"\
                         " for high-quality mate-pair library number <#> (<#> = 1,2,..,9)" + "\n")
        sys.stderr.write("--hqmp<#>-2\t<filename>\tfile with reverse reads"\
                         " for high-quality mate-pair library number <#> (<#> = 1,2,..,9)" + "\n")
        sys.stderr.write("--hqmp<#>-s\t<filename>\tfile with unpaired reads"\
                         " for high-quality mate-pair library number <#> (<#> = 1,2,..,9)" + "\n")
        sys.stderr.write("--hqmp<#>-<or>\torientation of reads"\
                         " for high-quality mate-pair library number <#> (<#> = 1,2,..,9; <or> = fr, rf, ff)" + "\n")
        sys.stderr.write("--nxmate<#>-1\t<filename>\tfile with forward reads"\
                             " for Lucigen NxMate library number <#> (<#> = 1,2,..,9)" + "\n")
        sys.stderr.write("--nxmate<#>-2\t<filename>\tfile with reverse reads"\
                             " for Lucigen NxMate library number <#> (<#> = 1,2,..,9)" + "\n")
        sys.stderr.write("--sanger\t<filename>\tfile with Sanger reads\n")
    if not mode == "rna":
        sys.stderr.write("--pacbio\t<filename>\tfile with PacBio reads\n")
        sys.stderr.write("--nanopore\t<filename>\tfile with Nanopore reads\n")
        sys.stderr.write("--tslr\t<filename>\tfile with TSLR-contigs\n")
    if not mode == "meta":
        sys.stderr.write("--trusted-contigs\t<filename>\tfile with trusted contigs\n")
        sys.stderr.write("--untrusted-contigs\t<filename>\tfile with untrusted contigs\n")

    if mode == "rna":
        sys.stderr.write("--ss-<type>\tstrand specific data, <type> = fr (normal) and rf (antisense)\n")

    sys.stderr.write("" + "\n")
    sys.stderr.write("Pipeline options:" + "\n")

    sys.stderr.write("--only-error-correction\truns only read error correction"\
                         " (without assembling)" + "\n")
    sys.stderr.write("--only-assembler\truns only assembling (without read error"\
                         " correction)" + "\n")

    if mode not in ["rna", "meta"]:
        sys.stderr.write("--careful\t\ttries to reduce number of mismatches and short indels" + "\n")
    sys.stderr.write("--continue\t\tcontinue run from the last available check-point" + "\n")
    if mode == "rna":
        sys.stderr.write("--restart-from\t<cp>\trestart run with updated options and from the specified check-point ('ec', 'as', 'last')" + "\n")
    else:
        sys.stderr.write("--restart-from\t<cp>\trestart run with updated options and from the specified check-point ('ec', 'as', 'k<int>', 'mc', 'last')" + "\n")
    sys.stderr.write("--disable-gzip-output\tforces error correction not to"\
                         " compress the corrected reads" + "\n")
    sys.stderr.write("--disable-rr\t\tdisables repeat resolution stage"\
                     " of assembling" + "\n")

    sys.stderr.write("" + "\n")
    sys.stderr.write("Advanced options:" + "\n")
    sys.stderr.write("--dataset\t<filename>\tfile with dataset description in YAML format" + "\n")
    if mode == "rna":
        sys.stderr.write("--fast\t\t\t\tspeeds up isoform detection, but may miss short and low-expressed isoforms\n")
    sys.stderr.write("-t/--threads\t<int>\t\tnumber of threads" + "\n")
    sys.stderr.write("\t\t\t\t[default: %s]\n" % THREADS)
    sys.stderr.write("-m/--memory\t<int>\t\tRAM limit for SPAdes in Gb"\
                         " (terminates if exceeded)" + "\n")
    sys.stderr.write("\t\t\t\t[default: %s]\n" % MEMORY)
    sys.stderr.write("--tmp-dir\t<dirname>\tdirectory for temporary files" + "\n")
    sys.stderr.write("\t\t\t\t[default: <output_dir>/tmp]" + "\n")
    if mode != 'rna':
        sys.stderr.write("-k\t\t<int,int,...>\tcomma-separated list of k-mer sizes" \
                         " (must be odd and" + "\n")
        sys.stderr.write("\t\t\t\tless than " + str(MAX_K + 1) + ") [default: 'auto']" + "\n")
    else:
        sys.stderr.write("-k\t\t<int>\t\tk-mer size (must be odd and less than " + str(MAX_K + 1) + ") " \
                         "[default: 'auto']\n")

    if mode not in ["rna", "meta"]:
        sys.stderr.write("--cov-cutoff\t<float>\t\tcoverage cutoff value (a positive float number, "
                         "or 'auto', or 'off') [default: 'off']" + "\n")
    sys.stderr.write("--phred-offset\t<33 or 64>\tPHRED quality offset in the"\
                         " input reads (33 or 64)" + "\n")
    sys.stderr.write("\t\t\t\t[default: auto-detect]" + "\n")    

    if show_hidden:
        sys.stderr.write("" + "\n")
        sys.stderr.write("HIDDEN options:" + "\n")
        sys.stderr.write("--debug\t\t\t\truns SPAdes in debug mode (keeps intermediate output)" + "\n")
        sys.stderr.write("--stop-after\t<cp>\truns SPAdes until the specified check-point ('ec', 'as', 'k<int>', 'mc') inclusive" + "\n")
        sys.stderr.write("--truseq\t\t\truns SPAdes in TruSeq mode\n")
        sys.stderr.write("--mismatch-correction\t\truns post processing correction"\
                             " of mismatches and short indels" + "\n")
        sys.stderr.write("--reference\t<filename>\tfile with reference for deep analysis"\
                             " (only in debug mode)" + "\n")
        sys.stderr.write("--series-analysis\t<filename>\tconfig for metagenomics-series-augmented reassembly" + "\n")
        sys.stderr.write("--configs-dir\t<configs_dir>\tdirectory with configs" + "\n")
        sys.stderr.write("-i/--iterations\t<int>\t\tnumber of iterations for read error"\
                             " correction [default: %s]\n" % ITERATIONS)
        sys.stderr.write("--read-buffer-size\t<int>\t\tsets size of read buffer for graph construction")
        sys.stderr.write("--bh-heap-check\t\t<value>\tsets HEAPCHECK environment variable"\
                             " for BayesHammer" + "\n")
        sys.stderr.write("--spades-heap-check\t<value>\tsets HEAPCHECK environment variable"\
                             " for SPAdes" + "\n")
        sys.stderr.write("--large-genome\tEnables optimizations for large genomes \n")
        sys.stderr.write("--save-gp\tEnables saving graph pack before repeat resolution (even without --debug) \n")
        sys.stderr.write("--hidden-cov-cutoff\t<float>\t\tcoverage cutoff value deeply integrated in simplification"\
                            " (a positive float number). Base coverage! Will be adjusted depending on K and RL! \n")
        sys.stderr.write("--read-cov-threshold\t<int>\t\tread median coverage threshold (non-negative integer)\n")
        sys.stderr.write("--help-hidden\tprints this usage message with all hidden options" + "\n")

    sys.stderr.flush()


def auto_K_allowed():
    return not k_mers and not single_cell and not iontorrent and not rna and not meta
    # kmers were set by default, not SC, not IonTorrent data and not rna and temporary not meta


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
    global cov_cutoff
    global tmp_dir
    global fast

    if threads is None:
        threads = THREADS
    if memory is None:
        if support.get_available_memory():
            memory = int(min(MEMORY, support.get_available_memory()))
        else:
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
    if cov_cutoff is None:
        cov_cutoff = 'off'
    if tmp_dir is None:
        tmp_dir = os.path.join(output_dir, TMP_DIR)
    if fast is None:
        fast = False


def set_test_options():
    global output_dir
    global single_cell
    global test_mode
    global meta

    output_dir = os.path.abspath('spades_test')
    single_cell = False
    meta = False
    test_mode = True


def save_restart_options(log):
    if dataset_yaml_filename:
        support.error("you cannot specify --dataset with --restart-from option!", log)
    if single_cell:
        support.error("you cannot specify --sc with --restart-from option!", log)
    if meta:
        support.error("you cannot specify --meta with --restart-from option!", log)
    if iontorrent:
        support.error("you cannot specify --iontorrent with --restart-from option!", log)
    if only_assembler:
        support.error("you cannot specify --only-assembler with --restart-from option!", log)
    if only_error_correction:
        support.error("you cannot specify --only-error-correction with --restart-from option!", log)
    if strand_specificity is not None:
        support.error("you cannot specify strand specificity (--ss-rf or --ss-fr) with --restart-from option!", log)

    global restart_k_mers
    global restart_careful
    global restart_mismatch_corrector
    global restart_disable_gzip_output
    global restart_disable_rr
    global restart_threads
    global restart_memory
    global restart_tmp_dir
    global restart_qvoffset
    global restart_cov_cutoff
    global restart_developer_mode
    global restart_reference
    global restart_configs_dir
    global restart_read_buffer_size
    global restart_fast

    restart_k_mers = k_mers
    restart_careful = careful
    restart_mismatch_corrector = mismatch_corrector
    restart_disable_gzip_output = disable_gzip_output
    restart_disable_rr = disable_rr
    restart_threads = threads
    restart_memory = memory
    restart_tmp_dir = tmp_dir
    restart_qvoffset = qvoffset
    restart_cov_cutoff = cov_cutoff
    restart_developer_mode = developer_mode
    restart_reference = reference
    restart_configs_dir = configs_dir
    restart_read_buffer_size = read_buffer_size
    restart_fast = fast


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
    global cov_cutoff
    global developer_mode
    global reference
    global configs_dir
    global read_buffer_size
    global original_k_mers
    global fast

    if restart_k_mers:
        original_k_mers = k_mers
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
    if restart_cov_cutoff is not None:
        cov_cutoff = restart_cov_cutoff
    if restart_developer_mode is not None:
        developer_mode = restart_developer_mode
    if restart_reference is not None:
        reference = restart_reference
    if restart_configs_dir is not None:
        configs_dir = restart_configs_dir
    if restart_read_buffer_size is not None:
        read_buffer_size = restart_read_buffer_size
    if restart_fast is not None:
        fast = restart_fast


def enable_truseq_mode():
    global truseq_mode
    global correct_scaffolds
    global run_truseq_postprocessing
    global K_MERS_SHORT
    global K_MERS_150
    global K_MERS_250
    global only_assembler
    global single_cell
    K_MERS_SHORT = [21,33,45,55]
    K_MERS_150 = [21,33,45,55,77]
    K_MERS_250 = [21,33,45,55,77,99,127]
    truseq_mode = True
    correct_scaffolds = True
    run_truseq_postprocessing = True
    only_assembler = True


def will_rerun(options):
    for opt, arg in options:
        if opt == '--continue' or opt.startswith('--restart-from'):  # checks both --restart-from k33 and --restart-from=k33
            return True
    return False
