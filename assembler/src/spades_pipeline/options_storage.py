import sys

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
dataset_yaml_filename = ''
threads = 16
memory = 250
tmp_dir = ''
k_mers = [21,33,55]
iterations = 1
qvoffset = None # auto-detect by default
developer_mode = False

# hidden options
mismatch_corrector = False
reference = ''
bh_heap_check = ''
spades_heap_check = ''
### END OF DEFAULT VALUES

# list of spades.py options
long_options = "12= threads= memory= tmp-dir= iterations= phred-offset= sc "\
               "only-error-correction only-assembler "\
               "disable-gzip-output help test debug reference= "\
               "bh-heap-check= spades-heap-check= help-hidden "\
               "config-file= dataset= mismatch-correction careful rectangles".split()
short_options = "o:1:2:s:k:t:m:i:h"

# adding multiple paired-end and mate-pair libraries support
reads_options = []
for i in xrange(MAX_LIBS_NUMBER):
    for type in ["pe", "mp"]:
        reads_options += ("%s%d-1= %s%d-2= %s%d-12= %s%d-s= %s%d-rf %s%d-fr %s%d-ff" % tuple([type, i + 1] * 7)).split()
long_options += reads_options
# for checking wether option corrseponds to reads or not
reads_options = map(lambda x:"--" + x.split('=')[0], reads_options)
reads_options += ["--12", "-1", "-2", "-s"]


def usage(spades_version, show_hidden=False):
    print >> sys.stderr, "SPAdes genome assembler v." + str(spades_version)
    print >> sys.stderr, "Usage:", sys.argv[0], "[options] -o <output_dir>"
    print >> sys.stderr, ""
    print >> sys.stderr, "Basic options:"
    print >> sys.stderr, "-o\t<output_dir>\tdirectory to store all the resulting files (required)"
    print >> sys.stderr, "--sc\t\t\tthis flag is required for MDA (single-cell)"\
                         " data"
    print >> sys.stderr, "--12\t<filename>\tfile with interlaced forward and reverse"\
                         " paired-end reads"
    print >> sys.stderr, "-1\t<filename>\tfile with forward paired-end reads"
    print >> sys.stderr, "-2\t<filename>\tfile with reverse paired-end reads"
    print >> sys.stderr, "-s\t<filename>\tfile with unpaired reads"
    print >> sys.stderr, "--pe<#>-12\t<filename>\tfile with interlaced"\
                         " reads for paired-end library number <#> (<#> = 1,2,3,4,5)"
    print >> sys.stderr, "--pe<#>-1\t<filename>\tfile with forward reads"\
                         " for paired-end library number <#> (<#> = 1,2,3,4,5)"
    print >> sys.stderr, "--pe<#>-2\t<filename>\tfile with reverse reads"\
                         " for paired-end library number <#> (<#> = 1,2,3,4,5)"
    print >> sys.stderr, "--pe<#>-s\t<filename>\tfile with unpaired reads"\
                         " for paired-end library number <#> (<#> = 1,2,3,4,5)"
    print >> sys.stderr, "--pe<#>-<or>\torientation of reads"\
                         " for paired-end library number <#> (<#> = 1,2,3,4,5; <or> = fr, rf, ff)"
    print >> sys.stderr, "--mp<#>-12\t<filename>\tfile with interlaced"\
                         " reads for mate-pair library number <#> (<#> = 1,2,3,4,5)"
    print >> sys.stderr, "--mp<#>-1\t<filename>\tfile with forward reads"\
                         " for mate-pair library number <#> (<#> = 1,2,3,4,5)"
    print >> sys.stderr, "--mp<#>-2\t<filename>\tfile with reverse reads"\
                         " for mate-pair library number <#> (<#> = 1,2,3,4,5)"
    print >> sys.stderr, "--mp<#>-s\t<filename>\tfile with unpaired reads"\
                         " for mate-pair library number <#> (<#> = 1,2,3,4,5)"
    print >> sys.stderr, "--mp<#>-<or>\torientation of reads"\
                         " for mate-pair library number <#> (<#> = 1,2,3,4,5; <or> = fr, rf, ff)"
    print >> sys.stderr, "--test\t\t\truns SPAdes on toy dataset"
    print >> sys.stderr, "-h/--help\t\tprints this usage message"

    print >> sys.stderr, ""
    print >> sys.stderr, "Pipeline options:"
    print >> sys.stderr, "--only-error-correction\truns only read error correction"\
                         " (without assembling)"
    print >> sys.stderr, "--only-assembler\truns only assembling (without read error"\
                         " correction)"
    print >> sys.stderr, "--disable-gzip-output\tforces error correction not to"\
                         " compress the corrected reads"
    print >> sys.stderr, "--careful\t\ttries to reduce number"\
                         " of mismatches and short indels"
    print >> sys.stderr, "--rectangles\t\tuses rectangle graph algorithm for repeat resolution"

    print >> sys.stderr, ""
    print >> sys.stderr, "Advanced options:"
    print >> sys.stderr, "--dataset\t<filename>\tfile with dataset description in YAML format"
    print >> sys.stderr, "-t/--threads\t<int>\t\tnumber of threads"
    print >> sys.stderr, "\t\t\t\t[default: %s]" % threads
    print >> sys.stderr, "-m/--memory\t<int>\t\tRAM limit for SPAdes in Gb"\
                         " (terminates if exceeded)"
    print >> sys.stderr, "\t\t\t\t[default: %s]" % memory
    print >> sys.stderr, "--tmp-dir\t<dirname>\tdirectory for read error correction"\
                         " temporary files"
    print >> sys.stderr, "\t\t\t\t[default: <output_dir>/corrected/tmp]"
    print >> sys.stderr, "-k\t\t<int,int,...>\tcomma-separated list of k-mer sizes"\
                         " (must be odd and"
    print >> sys.stderr, "\t\t\t\tless than 128) [default: " + ",".join(map(str, k_mers)) + "]"
    print >> sys.stderr, "-i/--iterations\t<int>\t\tnumber of iterations for read error"\
                         " correction"
    print >> sys.stderr, "\t\t\t\t[default: %s]" % iterations
    print >> sys.stderr, "--phred-offset\t<33 or 64>\tPHRED quality offset in the"\
                         " input reads (33 or 64)"
    print >> sys.stderr, "\t\t\t\t[default: auto-detect]"
    print >> sys.stderr, "--debug\t\t\t\truns SPAdes in debug mode (keeps intermediate output)"

    if show_hidden:
        print >> sys.stderr, ""
        print >> sys.stderr, "HIDDEN options:"
        print >> sys.stderr, "--mismatch-correction\t\truns post processing correction"\
                             " of mismatches and short indels"
        print >> sys.stderr, "--reference\t<filename>\tfile with reference for deep analysis"\
                             " (only in debug mode)"
        print >> sys.stderr, "--bh-heap-check\t\t<value>\tsets HEAPCHECK environment variable"\
                             " for BayesHammer"
        print >> sys.stderr, "--spades-heap-check\t<value>\tsets HEAPCHECK environment variable"\
                             " for SPAdes"
        print >> sys.stderr, "--help-hidden\tprints this usage message with all hidden options"


def set_test_options():
    global output_dir
    global single_cell

    output_dir = 'spades_test'
    single_cell = False