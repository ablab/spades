import sys

SUPPORTED_PYTHON_VERSIONS = ['2.4', '2.5', '2.6', '2.7', '3.2', '3.3']
# allowed reads extensions for BayesHammer and for thw whole SPAdes pipeline
BH_ALLOWED_READS_EXTENSIONS = ['.fq', '.fastq']
ALLOWED_READS_EXTENSIONS = BH_ALLOWED_READS_EXTENSIONS + ['.fa', '.fasta']
# reads could be gzipped
BH_ALLOWED_READS_EXTENSIONS += [x + '.gz' for x in BH_ALLOWED_READS_EXTENSIONS]
ALLOWED_READS_EXTENSIONS += [x + '.gz' for x in ALLOWED_READS_EXTENSIONS]

# we support up to MAX_LIBS_NUMBER paired-end libs and MAX_LIBS_NUMBER mate-pair libs
MAX_LIBS_NUMBER = 5

### DEFAULT VALUES:
# basic options
output_dir = ''
single_cell = False

# pipeline options
only_error_correction = False
only_assembler = False
disable_gzip_output = False
careful = False
rectangles = False

# advanced options
continue_mode = False
dataset_yaml_filename = ''
threads = 16
memory = 250
tmp_dir = ''
k_mers = None
k_mers_short = [21,33,55]
k_mers_150 = [21,33,55,77]
k_mers_250 = [21,33,55,77,99,127]
qvoffset = None # auto-detect by default
developer_mode = False

# hidden options
mismatch_corrector = False
reference = ''
iterations = 1
bh_heap_check = ''
spades_heap_check = ''
### END OF DEFAULT VALUES
dict_of_prefixes = dict()

# list of spades.py options
long_options = "12= threads= memory= tmp-dir= iterations= phred-offset= sc "\
               "only-error-correction only-assembler "\
               "disable-gzip-output help test debug reference= "\
               "bh-heap-check= spades-heap-check= help-hidden "\
               "config-file= dataset= mismatch-correction careful rectangles continue".split()
short_options = "o:1:2:s:k:t:m:i:h"

# adding multiple paired-end and mate-pair libraries support
reads_options = []
for i in range(MAX_LIBS_NUMBER):
    for type in ["pe", "mp"]:
        reads_options += ("%s%d-1= %s%d-2= %s%d-12= %s%d-s= %s%d-rf %s%d-fr %s%d-ff" % tuple([type, i + 1] * 7)).split()
long_options += reads_options
# for checking whether option corresponds to reads or not
reads_options = list(map(lambda x:"--" + x.split('=')[0], reads_options))
reads_options += ["--12", "-1", "-2", "-s"]


def usage(spades_version, show_hidden=False):
    sys.stderr.write("SPAdes genome assembler v." + str(spades_version) + "\n")
    sys.stderr.write("Usage: " + str(sys.argv[0]) + " [options] -o <output_dir>" + "\n")
    sys.stderr.write("" + "\n")
    sys.stderr.write("Basic options:" + "\n")
    sys.stderr.write("-o\t<output_dir>\tdirectory to store all the resulting files (required)" + "\n")
    sys.stderr.write("--sc\t\t\tthis flag is required for MDA (single-cell)"\
                         " data" + "\n")
    sys.stderr.write("--test\t\t\truns SPAdes on toy dataset" + "\n")
    sys.stderr.write("-h/--help\t\tprints this usage message" + "\n")

    sys.stderr.write("" + "\n")
    sys.stderr.write("Input data:" + "\n")
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

    sys.stderr.write("" + "\n")
    sys.stderr.write("Pipeline options:" + "\n")
    sys.stderr.write("--only-error-correction\truns only read error correction"\
                         " (without assembling)" + "\n")
    sys.stderr.write("--only-assembler\truns only assembling (without read error"\
                         " correction)" + "\n")
    sys.stderr.write("--careful\t\ttries to reduce number"\
                         " of mismatches and short indels" + "\n")
    sys.stderr.write("--continue\t\tcontinue run from the last available check-point" + "\n")
    sys.stderr.write("--disable-gzip-output\tforces error correction not to"\
                         " compress the corrected reads" + "\n")
    sys.stderr.write("--rectangles\t\tuses rectangle graph algorithm for repeat resolution" + "\n")

    sys.stderr.write("" + "\n")
    sys.stderr.write("Advanced options:" + "\n")
    sys.stderr.write("--dataset\t<filename>\tfile with dataset description in YAML format" + "\n")
    sys.stderr.write("-t/--threads\t<int>\t\tnumber of threads" + "\n")
    sys.stderr.write("\t\t\t\t[default: %s]\n" % threads)
    sys.stderr.write("-m/--memory\t<int>\t\tRAM limit for SPAdes in Gb"\
                         " (terminates if exceeded)" + "\n")
    sys.stderr.write("\t\t\t\t[default: %s]\n" % memory)
    sys.stderr.write("--tmp-dir\t<dirname>\tdirectory for read error correction"\
                         " temporary files" + "\n")
    sys.stderr.write("\t\t\t\t[default: <output_dir>/corrected/tmp]" + "\n")
    sys.stderr.write("-k\t\t<int,int,...>\tcomma-separated list of k-mer sizes"\
                         " (must be odd and" + "\n")
    sys.stderr.write("\t\t\t\tless than 128) [default: " + ",".join(map(str, k_mers_short)) + "]" + "\n")
    sys.stderr.write("--phred-offset\t<33 or 64>\tPHRED quality offset in the"\
                         " input reads (33 or 64)" + "\n")
    sys.stderr.write("\t\t\t\t[default: auto-detect]" + "\n")
    sys.stderr.write("--debug\t\t\t\truns SPAdes in debug mode (keeps intermediate output)" + "\n")

    if show_hidden:
        sys.stderr.write("" + "\n")
        sys.stderr.write("HIDDEN options:" + "\n")
        sys.stderr.write("--mismatch-correction\t\truns post processing correction"\
                             " of mismatches and short indels" + "\n")
        sys.stderr.write("--reference\t<filename>\tfile with reference for deep analysis"\
                             " (only in debug mode)" + "\n")
        sys.stderr.write("-i/--iterations\t<int>\t\tnumber of iterations for read error"\
                             " correction [default: %s]\n" % iterations)
        sys.stderr.write("--bh-heap-check\t\t<value>\tsets HEAPCHECK environment variable"\
                             " for BayesHammer" + "\n")
        sys.stderr.write("--spades-heap-check\t<value>\tsets HEAPCHECK environment variable"\
                             " for SPAdes" + "\n")
        sys.stderr.write("--help-hidden\tprints this usage message with all hidden options" + "\n")

    sys.stderr.flush()


def set_test_options():
    global output_dir
    global single_cell

    output_dir = 'spades_test'
    single_cell = False
