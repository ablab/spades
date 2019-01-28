/**
 * @file parasail_query
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * Reads fasta file of database sequences.
 * Optionally reads fasta file of query sequences.
 * Optionally filters inputs through suffix array.
 * Indexes input to learn length and end locations of each sequence.
 * Creates SA, LCP, and BWT.
 * Runs maximal pairs algorithm with the given minimum cutoff length.
 * For each sequence pair, performs alignment of user choice.
 * Output is csv of alignments between sequences.
 */
#include "config.h"

#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#if defined(HAVE_GETOPT) && defined(HAVE_UNISTD_H)
#include <unistd.h>
#elif defined(HAVE_WINDOWS_H)
#include "wingetopt/src/getopt.h"
#endif
#if defined(HAVE_POLL)
#include <err.h>
#include <poll.h>
#include <pwd.h>
#include <signal.h>
#endif
#if defined(HAVE_WINDOWS_H)
#include <windows.h>
#endif
#if defined(HAVE_FILELENGTH)
#include <io.h>
#endif

#include <algorithm>
#include <cctype>
#include <cfloat>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#define THREAD_DOC "      threads: system-specific default, must be >= 1\n"
#else
#define THREAD_DOC "      threads: Warning: ignored; OpenMP was not supported by your compiler\n"
#endif

#include "parasail.h"
#include "parasail/io.h"

#include "sais.h"

#if HAVE_VARIADIC_MACROS
#define eprintf(STREAM, ...) fprintf(STREAM, __VA_ARGS__); fflush(STREAM)
#else
#define eprintf fprintf
#endif

#define GB 1.0E-9

extern "C" size_t getMemorySize(void);

using ::std::bad_alloc;
using ::std::istringstream;
using ::std::make_pair;
using ::std::pair;
using ::std::set;
using ::std::size_t;
using ::std::stack;
using ::std::string;
using ::std::toupper;
using ::std::transform;
using ::std::vector;

typedef pair<int,int> Pair;

typedef set<Pair> PairSet;
typedef vector<Pair> PairVec;

struct quad {
    int lcp;
    int lb;
    int rb;
    vector<quad> children;

    quad()
        : lcp(0), lb(0), rb(INT_MAX), children() {}
    quad(int lcp, int lb, int rb)
        : lcp(lcp), lb(lb), rb(rb), children() {}
    quad(int lcp, int lb, int rb, vector<quad> children)
        : lcp(lcp), lb(lb), rb(rb), children(children) {}

    bool empty() { return rb == INT_MAX; }
};

inline static void pair_check(
        unsigned long &count_generated,
        PairSet &pairs,
        const int &i,
        const int &j,
        const int * const restrict SA,
        const unsigned char * const restrict BWT,
        const int * const restrict SID,
        const vector<int> &DB,
        const char &sentinal);

inline static void process(
        unsigned long &count_generated,
        PairSet &pairs,
        const quad &q,
        const int * const restrict SA,
        const unsigned char * const restrict BWT,
        const int * const restrict SID,
        const vector<int> &DB,
        const char &sentinal,
        const int &cutoff);

inline static void print_array(
        const char * filename_,
        const int * const restrict array,
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len);

inline static int self_score(
        const char * const restrict seq,
        int len,
        const parasail_matrix_t *matrix);

inline static void output_edges(
        FILE *fop,
        bool has_query,
        long sid_crossover,
        unsigned char *T,
        int AOL,
        int SIM,
        int OS,
        const parasail_matrix_t *matrix,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results,
        long long start,
        long long stop);

inline static void output_graph(
        FILE *fop,
        int which,
        unsigned char *T,
        int AOL,
        int SIM,
        int OS,
        const parasail_matrix_t *matrix,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results,
        vector<vector<pair<int,float> > > &graph,
        unsigned long &edge_count,
        long long start,
        long long stop);

inline static void output_stats(
        FILE *fop,
        bool has_query,
        long sid_crossover,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results,
        long long start,
        long long stop);

inline static void output_basic(
        FILE *fop,
        bool has_query,
        long sid_crossover,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results,
        long long start,
        long long stop);

inline static void output_emboss(
        FILE *fop,
        bool has_query,
        long sid_crossover,
        const parasail_matrix_t *matrix,
        const PairVec &vpairs,
        parasail_sequences_t *queries,
        parasail_sequences_t *sequences,
        const vector<parasail_result_t*> &results,
        long long start,
        long long stop);

inline static void output_ssw(
        FILE *fop,
        bool has_query,
        long sid_crossover,
        const parasail_matrix_t *matrix,
        const PairVec &vpairs,
        parasail_sequences_t *queries,
        parasail_sequences_t *sequences,
        const vector<parasail_result_t*> &results,
        long long start,
        long long stop);

inline static void output_sam(
        FILE *fop,
        bool use_sam_header,
        bool has_query,
        long sid_crossover,
        const parasail_matrix_t *matrix,
        const PairVec &vpairs,
        parasail_sequences_t *queries,
        parasail_sequences_t *sequences,
        const vector<parasail_result_t*> &results,
        long long start,
        long long stop);

inline static void output_trace(
        FILE *fop,
        bool has_query,
        long sid_crossover,
        unsigned char *T,
        const parasail_matrix_t *matrix,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results,
        long long start,
        long long stop);

inline static void output_tables(
        bool has_query,
        long sid_crossover,
        unsigned char *T,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results,
        long long start,
        long long stop);

inline static void output(
        bool is_stats,
        bool is_table,
        bool is_trace,
        bool edge_output,
        bool use_emboss_format,
        bool use_ssw_format,
        bool use_sam_format,
        bool use_sam_header,
        FILE *fop,
        bool has_query,
        long sid_crossover,
        unsigned char *T,
        int AOL,
        int SIM,
        int OS,
        const parasail_matrix_t *matrix,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        parasail_sequences_t *queries,
        parasail_sequences_t *sequences,
        const vector<parasail_result_t*> &results,
        long long start,
        long long stop);

inline static size_t parse_bytes(const char*);
static void set_signal_handler();

static void print_help(const char *progname, int status) {
    eprintf(stderr, "\nusage: %s "
            "[-a funcname] "
            "[-c cutoff] "
            "[-x] "
            "[-e gap_extend] "
            "[-o gap_open] "
            "[-m matrix] "
            "[-t threads] "
            "[-d] "
            "[-M match] "
            "[-X mismatch] "
            "[-k band size (for nw_banded)] "
            "[-l AOL] "
            "[-s SIM] "
            "[-i OS] "
            "[-v] "
            "[-V] "
            "-f file "
            "[-q query_file] "
            "[-g output_file] "
            "[-O output_format {EMBOSS,SAM,SAMH,SSW}] "
            "[-b batch_size] "
            "[-r memory_budget] "
            "\n\n",
            progname);
    eprintf(stderr, "Defaults:\n"
            "     funcname: sw_stats_striped_16\n"
            "       cutoff: 7, must be >= 1, exact match length cutoff\n"
            "           -x: if present, don't use suffix array filter\n"
            "   gap_extend: 1, must be >= 0\n"
            "     gap_open: 10, must be >= 0\n"
            "       matrix: blosum62\n"
            "           -d: if present, assume DNA alphabet\n"
            "        match: 1, must be >= 0\n"
            "     mismatch: 0, must be >= 0\n"
THREAD_DOC
            "          AOL: 80, must be 0 <= AOL <= 100, percent alignment length\n"
            "          SIM: 40, must be 0 <= SIM <= 100, percent exact matches\n"
            "           OS: 30, must be 0 <= OS <= 100, percent optimal score\n"
            "                                           over self score\n"
            "           -v: verbose output, report input parameters and timing\n"
            "           -V: verbose memory output, report memory use\n"
            "         file: no default, must be in FASTA format\n"
            "   query_file: no default, must be in FASTA format\n"
            "  output_file: parasail.csv\n"
            "output_format: no deafult, must be one of {EMBOSS,SAM,SAMH,SSW}\n"
            "   batch_size: 0 (calculate based on memory budget),\n"
            "               how many alignments before writing output\n"
            "memory_budget: 2GB or half available from system query (%.3f GB)\n",
        getMemorySize()/2.0*GB
            );
    exit(status);
}

#ifdef HAVE_POLL
static int stdin_has_data() {
    int timeout = 100; /* wait 100ms */
    struct pollfd fd;
    fd.fd = 0;
    fd.events = POLLIN;
    fd.revents = 0;
    int ret = poll(&fd, 1, timeout);
    return (ret > 0 && (fd.revents & POLLIN));
}
#elif defined(HAVE_WINDOWS_H)
#if !defined(HAVE_FILELENGTH)
#if defined(HAVE_STRUCT___STAT64) && defined(HAVE__FSTAT64)
#define STATBUF struct __stat64
#define FSTATFUNC _fstat64
#else
#define STATBUF struct stat
#define FSTATFUNC fstat
#endif
static long filelength(int fd) {
    int status;
    STATBUF buffer;
    status = FSTATFUNC(fd, &buffer);
    if (0 == status) {
        return buffer.st_size;
    }
    return -1;
}
#endif /* HAVE_FILELENGTH */
static int stdin_has_data() {
    int stdinHandle = fileno(stdin);
    long stdinFileLength = filelength(stdinHandle);
    if (stdinFileLength < 0) {
        return 0;
    }
    else if (stdinFileLength == 0) {
        int retval = fseek(stdin, 0, 0);
        if (retval < 0) {
            return 0;
        }
        else {
            return 1;
        }
    }
    return 1;
}
#endif

template <class info>
vector<long long> calc_batches(
        long long batch_size,
        bool verbose_memory,
        size_t memory_budget,
        const info *function_info,
        const PairVec &vpairs,
        const vector<long> &BEG,
        const vector<long> &END,
        size_t &memory_estimate)
{
    vector<long long> batches;
    batches.push_back(0);
    size_t result_size = sizeof(parasail_result_t);
    memory_estimate = 0;

    if (verbose_memory) {
        eprintf(stdout, "%20s: %.4f GB\n", "memory remaining", memory_budget*GB);
    }

    /* if user specified batch size on command line, we use it */
    if (0 != batch_size) {
        size_t batch_size_ = static_cast<size_t>(batch_size);
        for (size_t batch=batch_size_; batch<vpairs.size(); batch+=batch_size_) {
            batches.push_back(batch);
        }
        batches.push_back(vpairs.size());
        if (verbose_memory) {
            for (size_t batch=0; batch<batches.size()-1; ++batch) {
                long long start = batches[batch];
                long long stop = batches[batch+1];
                eprintf(stdout, "%20s: %lld-%lld\t%.4f GB\n", "batch",
                        start, stop, memory_estimate*GB);
            }
        }
        return batches;
    }

    if (function_info->is_table
            || function_info->is_rowcol
            || function_info->is_trace) {
        size_t multiplier = 1;
        if (function_info->is_table) {
            if (function_info->is_stats) {
                multiplier = 4;
                result_size += sizeof(parasail_result_extra_stats_t);
                result_size += sizeof(parasail_result_extra_stats_tables_t);
            }
            else {
                multiplier = 1;
                result_size += sizeof(parasail_result_extra_tables_t);
            }
        }
        else if (function_info->is_rowcol) {
            if (function_info->is_stats) {
                multiplier = 4;
                result_size += sizeof(parasail_result_extra_stats_t);
                result_size += sizeof(parasail_result_extra_stats_rowcols_t);
            }
            else {
                multiplier = 1;
                result_size += sizeof(parasail_result_extra_rowcols_t);
            }
        }
        else /* if (function_info->is_trace) */ {
            result_size += sizeof(parasail_result_extra_trace_t);
        }
        int i = vpairs[0].first;
        int j = vpairs[0].second;
        long i_beg = BEG[i];
        long i_end = END[i];
        long i_len = i_end-i_beg;
        long j_beg = BEG[j];
        long j_end = END[j];
        long j_len = j_end-j_beg;
        size_t current_size = 0;
        if (function_info->is_table) {
            current_size = sizeof(int) * multiplier * i_len * j_len;
        }
        else if (function_info->is_rowcol) {
            current_size = sizeof(int) * multiplier * (i_len + j_len);
        }
        else /* if (function_info->is_trace) */ {
            current_size = sizeof(int8_t) * multiplier * i_len * j_len;
        }
        memory_estimate = current_size;
        for (size_t index=1; index<vpairs.size(); ++index) {
            i = vpairs[index].first;
            j = vpairs[index].second;
            i_beg = BEG[i];
            i_end = END[i];
            i_len = i_end-i_beg;
            j_beg = BEG[j];
            j_end = END[j];
            j_len = j_end-j_beg;
            unsigned long local_size = 0;
            if (function_info->is_table) {
                local_size = sizeof(int) * multiplier * i_len * j_len;
            }
            else if (function_info->is_rowcol) {
                local_size = sizeof(int) * multiplier * (i_len + j_len);
            }
            else /* if (function_info->is_trace) */ {
                local_size = sizeof(int8_t) * multiplier * i_len * j_len;
            }
            if ((current_size + local_size) > memory_budget) {
                batches.push_back(index);
                if (current_size > memory_estimate) {
                    memory_estimate = current_size;
                }
                if (verbose_memory) {
                    long long start = batches[batches.size()-2];
                    long long stop = batches[batches.size()-1];
                    eprintf(stdout, "%20s: %lld-%lld\t%.4f GB\n", "batch",
                            start, stop, current_size*GB);
                }
                current_size = local_size;
            }
            else {
                current_size += local_size;
            }
        }
        if (current_size > memory_estimate) {
            memory_estimate = current_size;
        }
        batches.push_back(vpairs.size());
        if (verbose_memory) {
            long long start = batches[batches.size()-2];
            long long stop = batches[batches.size()-1];
            eprintf(stdout, "%20s: %lld-%lld\t%.4f GB\n", "batch",
                    start, stop, current_size*GB);
        }
    }
    else {
        if (function_info->is_stats) {
            result_size += sizeof(parasail_result_extra_stats_t);
        }
        /* the score, and perhaps stats */
        size_t results_per_batch = memory_budget / result_size;
        size_t how_many_batches = vpairs.size() / results_per_batch;
        for (size_t i=0; i<how_many_batches; ++i) {
            batches.push_back((i+1)*results_per_batch);
        }
        batches.push_back(vpairs.size());
        if (results_per_batch > vpairs.size()) {
            memory_estimate = result_size*vpairs.size();
        }
        else {
            memory_estimate = result_size*results_per_batch;
        }
        if (verbose_memory) {
            for (size_t batch=0; batch<batches.size()-1; ++batch) {
                long long start = batches[batch];
                long long stop = batches[batch+1];
                eprintf(stdout, "%20s: %lld-%lld\t%.4f GB\n", "batch",
                        start, stop, memory_estimate*GB);
            }
        }
    }

    return batches;
}


int main(int argc, char **argv) {
    FILE *fop = NULL;
    const char *fname = NULL;
    const char *qname = NULL;
    const char *oname = "parasail.csv";
    bool oname_from_user = false;
    parasail_sequences_t *sequences = NULL;
    parasail_sequences_t *queries = NULL;
    unsigned char *T = NULL;
    unsigned char *Q = NULL;
    int num_threads = -1;
    int *SA = NULL;
    int *LCP = NULL;
    unsigned char *BWT = NULL;
    int *SID = NULL;
    vector<long> BEG;
    vector<long> END;
    vector<int> DB;
    long n = 0;
    long t = 0;
    long q = 0;
    double start = 0;
    double finish = 0;
    long i = 0;
    long sid = 0;
    long sid_crossover = -1;
    char sentinal = 0;
    int cutoff = 7;
    bool use_filter = true;
    char *output_format = NULL;
    bool use_emboss_format = false;
    bool use_sam_format = false;
    bool use_sam_header = false;
    bool use_ssw_format = false;
    PairSet pairs;
    PairVec vpairs;
    unsigned long count_possible = 0;
    unsigned long count_generated = 0;
    unsigned long work = 0;
    int c = 0;
    const char *funcname = "sw_stats_striped_16";
    const parasail_function_info_t *function_info = NULL;
    const parasail_pfunction_info_t *pfunction_info = NULL;
    parasail_function_t *function = NULL;
    parasail_pfunction_t *pfunction = NULL;
    parasail_pcreator_t *pcreator = NULL;
    int is_banded = 0;
    int is_trace = 0;
    int kbandsize = 3;
    const char *matrixname = NULL;
    const parasail_matrix_t *matrix = NULL;
    int gap_open = 10;
    int gap_extend = 1;
    int match = 1;
    int mismatch = 0;
    bool use_dna = false;
    bool pairs_only = false;
    bool edge_output = false;
    bool graph_output = false;
    bool is_stats = true;
    bool is_table = false;
    bool has_query = false;
    const char *progname = "parasail_aligner";
    int AOL = 80;
    int SIM = 40;
    int OS = 30;
    bool verbose = false;
    bool verbose_memory = false;
    long long batch_size = 0;
    bool has_stdin = false;
    size_t memsize = 0; /* for temporary use */
    size_t memory_budget = 0;
    size_t bytes_used = 0;
    size_t profile_bits = 0;

    set_signal_handler();

    /* establish default memory budget */
    memory_budget = getMemorySize();
    if (0 == memory_budget) {
        /* error occurred, guess reasonable default */
        memory_budget = 2.0/GB;
    }
    else {
        /* use half available memory */
        memory_budget = memory_budget / 2;
    }

    /* Check arguments. */
    while ((c = getopt(argc, argv, "a:b:c:de:Ef:g:Ghi:k:l:m:M:o:O:pq:r:s:t:vVxX:")) != -1) {
        switch (c) {
            case 'a':
                funcname = optarg;
                break;
            case 'b':
                {
                    string numStr = optarg;
                    istringstream iss(numStr);
                    iss>>batch_size;
                    if (batch_size < 0) {
                    eprintf(stderr, "batch size must be >= 0\n");
                        print_help(progname, EXIT_FAILURE);
                    }
                }
                break;
            case 'c':
                cutoff = atoi(optarg);
                if (cutoff <= 0) {
                    eprintf(stderr, "cutoff must be > 0\n");
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 'd':
                use_dna = true;
                break;
            case 'e':
                gap_extend = atoi(optarg);
                if (gap_extend < 0) {
                    eprintf(stderr, "gap extend must be >= 0\n");
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 'E':
                edge_output = true;
                break;
            case 'f':
                fname = optarg;
                break;
            case 'g':
                oname = optarg;
                oname_from_user = true;
                break;
            case 'G':
                graph_output = true;
                break;
            case 'h':
                print_help(progname, EXIT_FAILURE);
                break;
            case 'i':
                OS = atoi(optarg);
                if (OS < 0 || OS > 100) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 'k':
                kbandsize = atoi(optarg);
                if (kbandsize <= 0) {
                    eprintf(stderr, "band size must be > 0\n");
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 'l':
                AOL = atoi(optarg);
                if (AOL < 0 || AOL > 100) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 'm':
                matrixname = optarg;
                break;
            case 'M':
                match = atoi(optarg);
                if (match < 0) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 'o':
                gap_open = atoi(optarg);
                if (gap_open < 0) {
                    eprintf(stderr, "gap open must be >= 0\n");
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 'O':
                output_format = optarg;
                break;
            case 'p':
                pairs_only = true;
                break;
            case 'q':
                qname = optarg;
                break;
            case 'r':
                memory_budget = parse_bytes(optarg);
                break;
            case 's':
                SIM = atoi(optarg);
                if (SIM < 0 || SIM > 100) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 't':
                num_threads = atoi(optarg);
#ifdef _OPENMP
#else
                eprintf(stdout, "-t number of threads requested, but OpenMP was not found during configuration. Running without threads.");
#endif
                break;
            case 'v':
                verbose = true;
                break;
            case 'V':
                verbose_memory = true;
                break;
            case 'x':
                use_filter = false;
                break;
            case 'X':
                mismatch = atoi(optarg);
                if (mismatch < 0) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case '?':
                if (optopt == 'a'
                        || optopt == 'c'
                        || optopt == 'e'
                        || optopt == 'f'
                        || optopt == 'g'
                        || optopt == 'i'
                        || optopt == 'k'
                        || optopt == 'l'
                        || optopt == 'm'
                        || optopt == 'M'
                        || optopt == 'o'
                        || optopt == 'q'
                        || optopt == 's'
                        || optopt == 't'
                        || optopt == 'X'
                        ) {
                    eprintf(stderr,
                            "Option -%c requires an argument.\n",
                            optopt);
                }
                else if (isprint(optopt)) {
                    eprintf(stderr, "Unknown option `-%c'.\n",
                            optopt);
                }
                else {
                    eprintf(stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                }
                print_help(progname, EXIT_FAILURE);
                break;
            default:
                eprintf(stderr, "default case in getopt\n");
                exit(EXIT_FAILURE);
        }
    }

    /* select the function */
    if (funcname) {
        if (NULL != strstr(funcname, "profile")) {
            pfunction_info = parasail_lookup_pfunction_info(funcname);
            if (NULL == pfunction_info) {
                eprintf(stderr, "Specified profile function not found.\n");
                exit(EXIT_FAILURE);
            }
            pfunction = pfunction_info->pointer;
            pcreator = pfunction_info->creator;
            if (NULL == strstr(funcname, "sat")) {
                profile_bits = atoi(pfunction_info->width);
            }
            else {
                profile_bits = 56; /* 8+16+32 */
            }
        }
        else {
            function_info = parasail_lookup_function_info(funcname);
            if (NULL == function_info && NULL != strstr(funcname, "nw_banded")) {
                is_banded = 1;
            }
            if (NULL == function_info && 0 == is_banded) {
                eprintf(stderr, "Specified function not found.\n");
                exit(EXIT_FAILURE);
            }
            function = function_info->pointer;
        }
    }
    else {
        eprintf(stderr, "No alignment function specified.\n");
        exit(EXIT_FAILURE);
    }

    has_stdin = stdin_has_data();
    if (has_stdin) {
        if (fname == NULL) {
            fname = "stdin";
        }
        else if (qname == NULL) {
            qname = "stdin";
            has_query = true;
        }
        else {
            eprintf(stderr, "input file, query file, and stdin detected; max inputs is 2\n");
            exit(EXIT_FAILURE);
        }
    }
    else {
        if (fname == NULL) {
            eprintf(stderr, "missing input file\n");
            print_help(progname, EXIT_FAILURE);
        }
        else if (qname != NULL) {
            has_query = true;
        }
    }

    is_trace = (NULL != strstr(funcname, "trace"));
    is_stats = (NULL != strstr(funcname, "stats"));
    is_table = (NULL != strstr(funcname, "table"));

    if (edge_output && graph_output) {
        eprintf(stderr, "Can only request one of edge or graph output.\n");
        exit(EXIT_FAILURE);
    }
    if ((edge_output || graph_output) && !is_stats) {
        eprintf(stderr, "Edge or graph output requested, but alignment function does not return statistics.\n");
        exit(EXIT_FAILURE);
    }
    if (graph_output && has_query) {
        eprintf(stderr, "Cannot specify a query file and output as a graph.\n");
        exit(EXIT_FAILURE);
    }

    if (NULL != output_format) {
        bool trace_warning = false;
        if (NULL != strstr(output_format, "SAMH")) {
            use_sam_format = true;
            use_sam_header = true;
            trace_warning = true;
        }
        else if (NULL != strstr(output_format, "SAM")) {
            use_sam_format = true;
            trace_warning = true;
        }
        else if (NULL != strstr(output_format, "EMBOSS")) {
            use_emboss_format = true;
            trace_warning = true;
        }
        else if (NULL != strstr(output_format, "SSW")) {
            use_ssw_format = true;
            trace_warning = true;
        }
        else {
            eprintf(stderr, "Unknown output format '%s'.\n", output_format);
            exit(EXIT_FAILURE);
        }
        if (!is_trace && trace_warning) {
            eprintf(stderr, "The selected output format '%s' requires an alignment function that returns a traceback.\n", output_format);
            exit(EXIT_FAILURE);
        }
    }
    else if (is_trace) {
        eprintf(stderr, "Please select trace output format.\n");
        exit(EXIT_FAILURE);
    }

    /* select the substitution matrix */
    if (NULL == matrixname && use_dna) {
        matrixname = "ACGT";
        matrix = parasail_matrix_create("ACGT", match, -mismatch);
    }
    else {
        if (NULL == matrixname) {
            matrixname = "blosum62";
        }
        matrix = parasail_matrix_lookup(matrixname);
        if (NULL == matrix) {
            /* try as a filename */
            matrix = parasail_matrix_from_file(matrixname);
        }
        if (NULL == matrix) {
            eprintf(stderr, "Specified substitution matrix not found.\n");
            exit(EXIT_FAILURE);
        }
    }

    if (is_trace && !oname_from_user) {
        oname = "stdout";
    }

    /* print the parameters for reference */
    if (verbose) {
        int major, minor, patch;
        parasail_version(&major, &minor, &patch);
        eprintf(stdout,
                "%20s: %d.%d.%d\n"
                "%20s: %s\n"
                "%20s: %d\n"
                "%20s: %s\n"
                "%20s: %d\n"
                "%20s: %d\n"
                "%20s: %s\n"
                "%20s: %d\n"
                "%20s: %d\n"
                "%20s: %d\n"
                "%20s: %s\n"
                "%20s: %s\n"
                "%20s: %s\n"
                "%20s: %lld\n"
                "%20s: %.4f GB\n",
                "parasail version", major, minor, patch,
                "funcname", funcname,
                "cutoff", cutoff,
                "use filter", use_filter ? "yes" : "no",
                "gap_extend", gap_extend,
                "gap_open", gap_open,
                "matrix", matrixname,
                "AOL", AOL,
                "SIM", SIM,
                "OS", OS,
                "file", fname,
                "query", (NULL == qname) ? "<no query>" : qname,
                "output", oname,
                "batch_size", batch_size,
                "memory_budget", memory_budget*GB
                    );
        if (use_dna) {
            eprintf(stdout,
                    "%20s: %d\n"
                    "%20s: %d\n",
                    "match", match,
                    "mismatch", mismatch);
        }
    }
    if (verbose_memory) {
        eprintf(stdout, "%20s: %.4f GB\n", "memory budget", memory_budget*GB);
    }

    /* Best to know early whether we can open the output file. */
    if (oname_from_user || !is_trace) {
        if ((fop = fopen(oname, "w")) == NULL) {
            eprintf(stderr, "%s: Cannot open output file `%s': ",
                    progname, oname);
            perror("fopen");
            exit(EXIT_FAILURE);
        }
    }
    else {
        fop = stdout;
    }

    start = parasail_time();
    if (!has_query) {
        size_t count = 0;
        sequences = parasail_sequences_from_file(fname);
        T = (unsigned char*)parasail_sequences_pack(sequences, &count);
        n = count;
        if (is_trace) {
            /* This does not include name, comment, or qual strings. */
            bytes_used += sizeof(parasail_sequence_t) * sequences->l;
            bytes_used += sequences->characters;
        }
        else {
            /* If we aren't using tracebacks, we can discard sequences now. */
            parasail_sequences_free(sequences);
        }
        bytes_used += count;
    }
    else {
        size_t count = 0;
        sequences = parasail_sequences_from_file(fname);
        T = (unsigned char*)parasail_sequences_pack(sequences, &count);
        t = count;
        if (is_trace) {
            /* This does not include name, comment, or qual strings. */
            bytes_used += sizeof(parasail_sequence_t) * sequences->l;
            bytes_used += sequences->characters;
        }
        else {
            /* If we aren't using tracebacks, we can discard sequences now. */
            parasail_sequences_free(sequences);
        }
        bytes_used += count;
        queries = parasail_sequences_from_file(qname);
        Q = (unsigned char*)parasail_sequences_pack(queries, &count);
        q = count;
        if (is_trace) {
            /* This does not include name, comment, or qual strings. */
            bytes_used += sizeof(parasail_sequence_t) * queries->l;
            bytes_used += queries->characters;
        }
        else {
            /* If we aren't using tracebacks, we can discard queries now. */
            parasail_sequences_free(queries);
        }
        bytes_used += count;
        n = t+q;
        /* realloc T and copy Q into it */
        T = (unsigned char*)realloc(T, (n+1)*sizeof(unsigned char));
        if (T == NULL) {
            eprintf(stderr, "%s: Cannot reallocate memory.\n", progname);
            perror("realloc");
            exit(EXIT_FAILURE);
        }
        (void)memcpy(T+t, Q, q);
        free(Q);
    }
    T[n] = '\0';
    finish = parasail_time();
    if (verbose) {
        eprintf(stdout, "%20s: %.4f seconds\n", "read and pack time", finish-start);
    }
    if (verbose_memory) {
        eprintf(stdout, "%20s: %.4f GB\n", "read and pack memory", bytes_used*GB);
    }

    /* Allocate memory for sequence ID array. */
    if (use_filter) {
        memsize = (size_t)n * sizeof(int);
        SID = (int *)malloc(memsize);
        if(SID == NULL) {
            perror("malloc");
            eprintf(stderr, "%s: Cannot allocate suffix ID memory.\n", progname);
            eprintf(stderr, "Attempted %llu bytes. %llu already used.\n",
                    (unsigned long long)memsize,
                    (unsigned long long)bytes_used);
            exit(EXIT_FAILURE);
        }
        bytes_used += memsize;
    }

    /* determine sentinal */
    if (sentinal == 0) {
        long off = 0;
        while (!isgraph(T[n-off])) {
            ++off;
        }
        sentinal = T[n-off];
    }
    if (verbose) {
        eprintf(stdout, "%20s: %c\n", "sentinal", sentinal);
    }

    /* determine actual end of file (last char) */
    {
        long off = 0;
        while (!isgraph(T[n-off])) {
            ++off;
        }
        n = n - off + 1;
    }
    if (verbose) {
        eprintf(stdout, "%20s: %ld\n", "end of packed buffer", n);
    }

    /* scan T from left to count number of sequences */
    sid = 0;
    for (i=0; i<n; ++i) {
        if (T[i] == sentinal) {
            ++sid;
        }
    }
    if (0 == sid) { /* no sentinal found */
        eprintf(stderr, "no sentinal(%c) found in input\n", sentinal);
        exit(EXIT_FAILURE);
    }
    if (verbose) {
        eprintf(stdout, "%20s: %ld\n", "number of sequences", sid);
    }

    /* scan T from left to build sequence ID and end index */
    /* allocate vectors now that number of sequences is known */
    try {
        memsize = sizeof(long)*(sid+1);
        BEG.reserve(sid+1);
        bytes_used += memsize;
        END.reserve(sid+1);
        bytes_used += memsize;
        if (use_filter) {
            memsize = sizeof(int)*(sid+1);
            DB.reserve(sid+1);
            bytes_used += memsize;
        }
    } catch (const bad_alloc&) {
        eprintf(stderr, "Cannot allocate memory for vectors.\n");
        eprintf(stderr, "Attempted %llu bytes. %llu already used.\n",
                (unsigned long long)memsize,
                (unsigned long long)bytes_used);
        exit(EXIT_FAILURE);
    }
    sid = 0;
    BEG.push_back(0);
    if (use_filter) {
        /* look for mixed case also */
        char found_lower = 0;
        char found_upper = 0;
        for (i=0; i<n; ++i) {
            SID[i] = sid;
            if (T[i] == sentinal) {
                END.push_back(i);
                BEG.push_back(i+1);
                DB.push_back(i<t);
                if (-1 == sid_crossover && i>=t) {
                    sid_crossover = sid;
                }
                ++sid;
            }
            found_lower = found_lower || (T[i] >= 'a' && T[i] <= 'z');
            found_upper = found_upper || (T[i] >= 'A' && T[i] <= 'Z');
        }
        if (found_lower && found_upper) {
            /* make the whole thing upper case */
            for (i=0; i<n; ++i) {
                T[i] = toupper(T[i]);
            }
        }
    }
    else {
        for (i=0; i<n; ++i) {
            if (T[i] == sentinal) {
                END.push_back(i);
                BEG.push_back(i+1);
                if (-1 == sid_crossover && i>=t) {
                    sid_crossover = sid;
                }
                ++sid;
            }
        }
    }

    /* if we don't have a query file, clear the DB flags */
    if (!has_query) {
        /* vector::clear() might not deallocate memory.
         * Use swap with temporary instead */
        vector<int>().swap(DB);
        sid_crossover = -1;
    }
    else if (verbose) {
        eprintf(stdout, "%20s: %ld\n", "number of queries", sid-sid_crossover);
        eprintf(stdout, "%20s: %ld\n", "number of db seqs", sid_crossover);
    }

    /* use the enhanced SA filter */
    if (use_filter) {
        size_t memsize_local = 0;
        /* Allocate memory for enhanced SA. */
        memsize = (size_t)(n+1) * sizeof(int); /* +1 for LCP */
        SA = (int *)malloc(memsize);
        if (SA == NULL) {
            perror("malloc");
            eprintf(stderr, "%s: Cannot allocate SA memory.\n", progname);
            eprintf(stderr, "Attempted %llu bytes. %llu already used.\n",
                    (unsigned long long)memsize,
                    (unsigned long long)bytes_used);
            exit(EXIT_FAILURE);
        }
        memsize_local += memsize;
        bytes_used += memsize;
        memsize = (size_t)(n+1) * sizeof(int); /* +1 for lcp tree */
        LCP = (int *)malloc(memsize);
        if (LCP == NULL) {
            perror("malloc");
            eprintf(stderr, "%s: Cannot allocate LCP memory.\n", progname);
            eprintf(stderr, "Attempted %llu bytes. %llu already used.\n",
                    (unsigned long long)memsize,
                    (unsigned long long)bytes_used);
            exit(EXIT_FAILURE);
        }
        memsize_local += memsize;
        bytes_used += memsize;
        memsize = (size_t)(n+1) * sizeof(unsigned char);
        BWT = (unsigned char *)malloc(memsize);
        if (BWT == NULL) {
            perror("malloc");
            eprintf(stderr, "%s: Cannot allocate BWT memory.\n", progname);
            eprintf(stderr, "Attempted %llu bytes. %llu already used.\n",
                    (unsigned long long)memsize,
                    (unsigned long long)bytes_used);
            exit(EXIT_FAILURE);
        }
        memsize_local += memsize;
        bytes_used += memsize;

        /* Construct the suffix and LCP arrays.
         * The following sais routine is from Fischer, with bugs fixed. */
        start = parasail_time();
        if(sais(T, SA, LCP, (int)n) != 0) {
            eprintf(stderr, "%s: Cannot allocate memory.\n", progname);
            exit(EXIT_FAILURE);
        }
        finish = parasail_time();
        if (verbose) {
            eprintf(stdout,"%20s: %.4f seconds\n", "induced SA time", finish-start);
        }

        /* construct naive BWT: */
        start = parasail_time();
        for (i = 0; i < n; ++i) {
            BWT[i] = (SA[i] > 0) ? T[SA[i]-1] : sentinal;
        }
        finish = parasail_time();
        if (verbose) {
            eprintf(stdout, "%20s: %.4f seconds\n", "naive BWT time", finish-start);
        }

        /* "fix" the LCP array to clamp LCP's that are too long */
        start = parasail_time();
        for (i = 0; i < n; ++i) {
            int len = END[SID[SA[i]]] - SA[i]; /* don't include sentinal */
            if (LCP[i] > len) LCP[i] = len;
        }
        finish = parasail_time();
        if (verbose) {
            eprintf(stdout, "%20s: %.4f seconds\n", "clamp LCP time", finish-start);
        }

        /* The GSA we create will put all sentinals either at the beginning
         * or end of the SA. We don't want to count all of the terminals,
         * nor do we want to process them in our bottom-up traversal. */
        /* do the sentinals appear at the beginning or end of SA? */
        int bup_start = 1;
        int bup_stop = n;
        if (T[SA[0]] == sentinal) {
            /* sentinals at beginning */
            bup_start = sid+1;
            bup_stop = n;
        }
        else if (T[SA[n-1]] == sentinal) {
            /* sentinals at end */
            bup_start = 1;
            bup_stop = n-sid;
        }
        else {
            eprintf(stderr, "sentinals not found at beginning or end of SA\n");
            exit(EXIT_FAILURE);
        }

        /* DFS of enhanced SA, from Abouelhoda et al */
        start = parasail_time();
        count_generated = 0;
        LCP[n] = 0; /* doesn't really exist, but for the root */
        {
            stack<quad> the_stack;
            quad last_interval;
            the_stack.push(quad());
            for (i = bup_start; i <= bup_stop; ++i) {
                int lb = i - 1;
                while (LCP[i] < the_stack.top().lcp) {
                    the_stack.top().rb = i - 1;
                    last_interval = the_stack.top();
                    the_stack.pop();
                    process(count_generated, pairs, last_interval, SA, BWT, SID, DB, sentinal, cutoff);
                    lb = last_interval.lb;
                    if (LCP[i] <= the_stack.top().lcp) {
                        last_interval.children.clear();
                        the_stack.top().children.push_back(last_interval);
                        last_interval = quad();
                    }
                }
                if (LCP[i] > the_stack.top().lcp) {
                    if (!last_interval.empty()) {
                        last_interval.children.clear();
                        the_stack.push(quad(LCP[i],lb,INT_MAX,vector<quad>(1, last_interval)));
                        last_interval = quad();
                    }
                    else {
                        the_stack.push(quad(LCP[i],lb,INT_MAX));
                    }
                }
            }
            the_stack.top().rb = bup_stop - 1;
            process(count_generated, pairs, the_stack.top(), SA, BWT, SID, DB, sentinal, cutoff);
        }
        finish = parasail_time();
        if (!has_query) {
            count_possible = ((unsigned long)sid)*((unsigned long)sid-1)/2;
        } else {
            count_possible = (sid-sid_crossover)*sid_crossover;
        }
        if (verbose) {
            eprintf(stdout, "%20s: %.4f seconds\n", "ESA time", finish-start);
            eprintf(stdout, "%20s: %lu\n", "possible pairs", count_possible);
            eprintf(stdout, "%20s: %lu\n", "generated pairs", count_generated);
        }

        /* Deallocate memory. */
        free(SID);
        free(SA);
        free(LCP);
        free(BWT);
        bytes_used -= (size_t)n * sizeof(int); /* SID */
        bytes_used -= memsize_local; /* SA,LCP,BWT */

        if (verbose) {
            eprintf(stdout, "%20s: %zu\n", "unique pairs", pairs.size());
        }
    }
    else {
        /* don't use enhanced SA filter -- generate all pairs */
        start = parasail_time();
        if (!has_query) {
            /* no query file, so all against all comparison */
            for (int i=0; i<sid; ++i) {
                for (int j=i+1; j<sid; ++j) {
                    vpairs.push_back(make_pair(i,j));
                }
            }
        }
        else {
            /* query given, so only compare query against database */
            for (int i=sid_crossover; i<sid; ++i) {
                for (int j=0; j<sid_crossover; ++j) {
                    vpairs.push_back(make_pair(i,j));
                }
            }
        }
        finish = parasail_time();
        if (verbose) {
            eprintf(stdout, "%20s: %.4f seconds\n", "enumerate time", finish-start);
            eprintf(stdout, "%20s: %zu\n", "unique pairs", vpairs.size());
        }
    }

    if (!is_trace && pairs_only) {
        /* Done with input text. */
        free(T);
        if (vpairs.empty() && !pairs.empty()) {
            for (PairSet::iterator it=pairs.begin(); it!=pairs.end(); ++it) {
                int i = it->first;
                int j = it->second;
                eprintf(fop, "%d,%d\n", i, j);
            }
        }
        else if (!vpairs.empty() && pairs.empty()) {
            for (PairVec::iterator it=vpairs.begin(); it!=vpairs.end(); ++it) {
                int i = it->first;
                int j = it->second;
                eprintf(fop, "%d,%d\n", i, j);
            }
        }
        else {
            eprintf(stderr, "pairs and vpairs were empty\n");
            exit(EXIT_FAILURE);
        }
        fclose(fop);
        return 0;
    }

#ifdef _OPENMP
    if (-1 == num_threads) {
        num_threads = omp_get_max_threads();
    }
    else if (num_threads >= 1) {
        omp_set_num_threads(num_threads);
    }
    else {
        eprintf(stderr, "invalid number of threads chosen (%d)\n", num_threads);
        exit(EXIT_FAILURE);
    }
    if (verbose) {
        eprintf(stdout, "%20s: %d\n", "omp num threads", num_threads);
    }
#endif

    /* OpenMP can't iterate over an STL set. Convert to STL vector. */
    start = parasail_time();
    if (vpairs.empty()) {
        if (pairs.empty()) {
            if (use_filter) {
                eprintf(stderr, "no alignment work, either the filter removed all alignemnts or the input file(s) were empty\n");
                exit(EXIT_FAILURE);
            }
            else {
                eprintf(stderr, "no alignment work, perhaps the input file(s) were empty\n");
                exit(EXIT_FAILURE);
            }
        }
        vpairs.assign(pairs.begin(), pairs.end());
        /* set::clear() might not deallocate memory.
         * Use swap with temporary instead */
        PairSet().swap(pairs);
        /* don't increment bytes_used since memory use should remain same */
    }
    if (!pairs.empty()) {
        eprintf(stderr, "failed to free pair memory, continuing\n");
    }
    /* finally tally the pair memory */
    bytes_used += vpairs.size()*sizeof(Pair);
    /* pre-allocate result pointers */
    vector<parasail_result_t*> results(vpairs.size(),
            static_cast<parasail_result_t*>(NULL));
    bytes_used += vpairs.size()*sizeof(parasail_result_t*);
    finish = parasail_time();
    if (verbose) {
        eprintf(stdout, "%20s: %.4f seconds\n", "openmp prep time", finish-start);
    }
    if (verbose_memory) {
        eprintf(stdout, "%20s: %.4f GB\n", "openmp prep memory", bytes_used*GB);
    }

    /* create profiles, if necessary */
    vector<parasail_profile_t*> profiles;
    if (pfunction) {
        start = parasail_time();
        set<int> profile_indices_set;
        for (size_t index=0; index<vpairs.size(); ++index) {
            profile_indices_set.insert(vpairs[index].first);
        }
        vector<int> profile_indices(
                profile_indices_set.begin(),
                profile_indices_set.end());
        profiles.assign(sid, static_cast<parasail_profile_t*>(NULL));
        bytes_used += sizeof(parasail_profile_t*) * sid;
        finish = parasail_time();
        if (verbose) {
            eprintf(stdout, "%20s: %.4f seconds\n", "profile init", finish-start);
        }
        start = parasail_time();
#pragma omp parallel for schedule(guided)
        for (long long index=0; index<(long long)profile_indices.size(); ++index)
        {
            int i = profile_indices[index];
            long i_beg = BEG[i];
            long i_end = END[i];
            long i_len = i_end-i_beg;
            size_t local_mem = matrix->size * i_len * profile_bits;
            profiles[i] = pcreator((const char*)&T[i_beg], i_len, matrix);
#pragma omp atomic
            bytes_used += local_mem;
        }
        finish = parasail_time();
        if (verbose) {
            eprintf(stdout, "%20s: %.4f seconds\n", "profile creation", finish-start);
        }
        if (verbose_memory) {
            eprintf(stdout, "%20s: %.4f GB\n", "profile creation memory", bytes_used*GB);
        }
    }

    if (bytes_used > memory_budget) {
        eprintf(stderr, "memory budget exceeded prior to alignment phase\n");
        return 0;
    }

    /* align pairs */
    start = parasail_time();
    if (function) {
        size_t memory_estimate;
        long long vpairs_size = (long long)vpairs.size();
        vector<long long> batches = calc_batches(
                batch_size,
                verbose && verbose_memory,
                memory_budget-bytes_used, function_info,
                vpairs, BEG, END,
                memory_estimate);
        bytes_used += memory_estimate;
        vector<vector<pair<int,float> > > graph;
        unsigned long edge_count = 0;
        if (graph_output) {
            graph.resize(sid);
        }
        for (size_t batch=0; batch<batches.size()-1; ++batch) {
            long long start = batches[batch];
            long long stop = batches[batch+1];
            if (stop > vpairs_size) stop = vpairs_size;
#pragma omp parallel for schedule(guided)
            for (long long index=start; index<stop; ++index)
            {
                int i = vpairs[index].first;
                int j = vpairs[index].second;
                long i_beg = BEG[i];
                long i_end = END[i];
                long i_len = i_end-i_beg;
                long j_beg = BEG[j];
                long j_end = END[j];
                long j_len = j_end-j_beg;
                unsigned long local_work = i_len * j_len;
                parasail_result_t *result = function(
                        (const char*)&T[i_beg], i_len,
                        (const char*)&T[j_beg], j_len,
                        gap_open, gap_extend, matrix);
#pragma omp atomic
                work += local_work;
                results[index] = result;
            }
            if (graph_output) {
                output_graph(NULL, 0, T, AOL, SIM, OS, matrix, BEG,
                        END, vpairs, results, graph, edge_count, start,
                        stop);
            }
            else {
                output(is_stats, is_table, is_trace, edge_output,
                        use_emboss_format, use_ssw_format,
                        use_sam_format, use_sam_header, fop, has_query,
                        sid_crossover, T, AOL, SIM, OS, matrix,
                        BEG, END, vpairs, queries, sequences, results,
                        start, stop);
            }
            for (long long index=start; index<stop; ++index) {
                parasail_result_t *result = results[index];
                parasail_result_free(result);
            }
        }
        if (graph_output) {
            output_graph(fop, 0, T, AOL, SIM, OS, matrix, BEG,
                    END, vpairs, results, graph, edge_count, 0,
                    vpairs_size);
        }
    }
    else if (is_banded) {
        size_t memory_estimate;
        long long vpairs_size = (long long)vpairs.size();
        vector<long long> batches = calc_batches(
                batch_size,
                verbose && verbose_memory,
                memory_budget-bytes_used, function_info,
                vpairs, BEG, END,
                memory_estimate);
        bytes_used += memory_estimate;
        for (size_t batch=0; batch<batches.size()-1; ++batch) {
            long long start = batches[batch];
            long long stop = batches[batch+1];
            if (stop > vpairs_size) stop = vpairs_size;
#pragma omp parallel for schedule(guided)
            for (long long index=start; index<stop; ++index)
            {
                int i = vpairs[index].first;
                int j = vpairs[index].second;
                long i_beg = BEG[i];
                long i_end = END[i];
                long i_len = i_end-i_beg;
                long j_beg = BEG[j];
                long j_end = END[j];
                long j_len = j_end-j_beg;
                unsigned long local_work = i_len * j_len;
                parasail_result_t *result = parasail_nw_banded(
                        (const char*)&T[i_beg], i_len,
                        (const char*)&T[j_beg], j_len,
                        gap_open, gap_extend, kbandsize, matrix);
#pragma omp atomic
                work += local_work;
                results[index] = result;
            }
            output(is_stats, is_table, is_trace, edge_output,
                    use_emboss_format, use_ssw_format, use_sam_format,
                    use_sam_header, fop, has_query, sid_crossover, T, AOL,
                    SIM, OS, matrix, BEG, END, vpairs, queries, sequences,
                    results, start, stop);
            for (long long index=start; index<stop; ++index) {
                parasail_result_t *result = results[index];
                parasail_result_free(result);
            }
        }
    }
    else if (pfunction) {
        size_t memory_estimate;
        long long vpairs_size = (long long)vpairs.size();
        vector<long long> batches = calc_batches(
                batch_size,
                verbose && verbose_memory,
                memory_budget-bytes_used, pfunction_info,
                vpairs, BEG, END,
                memory_estimate);
        bytes_used += memory_estimate;
        vector<vector<pair<int,float> > > graph;
        unsigned long edge_count = 0;
        if (graph_output) {
            graph.resize(sid);
        }
        for (size_t batch=0; batch<batches.size()-1; ++batch) {
            long long start = batches[batch];
            long long stop = batches[batch+1];
            if (stop > vpairs_size) stop = vpairs_size;
#pragma omp parallel for schedule(guided)
            for (long long index=start; index<stop; ++index)
            {
                int i = vpairs[index].first;
                int j = vpairs[index].second;
                long j_beg = BEG[j];
                long j_end = END[j];
                long j_len = j_end-j_beg;
                parasail_profile_t *profile = profiles[i];
                if (NULL == profile) {
                    eprintf(stderr, "BAD PROFILE %d\n", i);
                    exit(EXIT_FAILURE);
                }
                unsigned long local_work = profile->s1Len * j_len;
                parasail_result_t *result = pfunction(
                        profile, (const char*)&T[j_beg], j_len,
                        gap_open, gap_extend);
#pragma omp atomic
                work += local_work;
                results[index] = result;
            }
            if (graph_output) {
                output_graph(NULL, 0, T, AOL, SIM, OS, matrix, BEG,
                        END, vpairs, results, graph, edge_count, start,
                        stop);
            }
            else {
                output(is_stats, is_table, is_trace, edge_output,
                        use_emboss_format, use_ssw_format,
                        use_sam_format, use_sam_header, fop, has_query,
                        sid_crossover, T, AOL, SIM, OS, matrix,
                        BEG, END, vpairs, queries, sequences, results,
                        start, stop);
            }
            for (long long index=start; index<stop; ++index) {
                parasail_result_t *result = results[index];
                parasail_result_free(result);
            }
        }
        if (graph_output) {
            output_graph(fop, 0, T, AOL, SIM, OS, matrix, BEG,
                    END, vpairs, results, graph, edge_count, 0,
                    vpairs_size);
        }
    }
    else {
        /* shouldn't get here */
        eprintf(stderr, "alignment function was not properly set (shouldn't happen)\n");
        exit(EXIT_FAILURE);
    }
    finish = parasail_time();
    if (verbose) {
        eprintf(stdout, "%20s: %lu cells\n", "work", work);
        eprintf(stdout, "%20s: %.4f seconds\n", "alignment time", finish-start);
        eprintf(stdout, "%20s: %.4f \n", "gcups", double(work)/(finish-start)/1000000000);
    }
    if (verbose_memory) {
        eprintf(stdout, "%20s: %.4f GB\n", "post-result memory", bytes_used*GB);
    }

    if (pfunction) {
        start = parasail_time();
#pragma omp parallel for schedule(guided)
        for (long long index=0; index<(long long)profiles.size(); ++index)
        {
            if (NULL != profiles[index]) {
                parasail_profile_free(profiles[index]);
            }
        }
        profiles.clear();
        finish = parasail_time();
        if (verbose) {
            eprintf(stdout, "%20s: %.4f seconds\n", "profile cleanup", finish-start);
        }
    }

    /* close output file */
    if (oname_from_user || !is_trace) {
        fclose(fop);
    }

    /* Done with input text. */
    free(T);

    /* Done with sequences if we were using tracebacks. */
    if (is_trace) {
        parasail_sequences_free(sequences);
        if (has_query) {
            parasail_sequences_free(queries);
        }
    }

    return 0;
}


inline static void pair_check(
        unsigned long &count_generated,
        PairSet &pairs,
        const int &i,
        const int &j,
        const int * const restrict SA,
        const unsigned char * const restrict BWT,
        const int * const restrict SID,
        const vector<int> &DB,
        const char &sentinal)
{
    const int &sidi = SID[SA[i]];
    const int &sidj = SID[SA[j]];
    if (BWT[i] != BWT[j] || BWT[i] == sentinal) {
        if (DB.empty()) {
            if (sidi != sidj) {
                ++count_generated;
                if (sidi < sidj) {
                    pairs.insert(make_pair(sidi,sidj));
                }
                else {
                    pairs.insert(make_pair(sidj,sidi));
                }
            }
        }
        else {
            if (sidi != sidj && DB[sidi] != DB[sidj]) {
                ++count_generated;
                if (sidi > sidj) {
                    pairs.insert(make_pair(sidi,sidj));
                }
                else {
                    pairs.insert(make_pair(sidj,sidi));
                }
            }
        }
    }
}

/* try to reduce number of duplicate pairs generated */
/* we observe that l-intervals (i.e. internal nodes) always have at
 * least two children, but these children could be singleton
 * l-intervals, e.g., [i..j]=[1..1], in addition to l-intervals with
 * non-singleton ranges/quads. For each l-interval, we take the cross
 * product of its child l-intervals. Naively, we could take the cross
 * product of the entire lb/rb range of the l-interval, but this
 * generates too many duplicate pairs. Instead, the complexity should be
 * bounded by the number of exact matches...
 */
inline static void process(
        unsigned long &count_generated,
        PairSet &pairs,
        const quad &q,
        const int * const restrict SA,
        const unsigned char * const restrict BWT,
        const int * const restrict SID,
        const vector<int> &DB,
        const char &sentinal,
        const int &cutoff)
{
    const size_t n_children = q.children.size();
    size_t child_index = 0;

    if (q.lcp < cutoff) return;

    if (n_children) {
        for (int i=q.lb; i<=q.rb; ++i) {
            int j = i+1;
            if (child_index < n_children) {
                if (i >= q.children[child_index].lb) {
                    j = q.children[child_index].rb+1;
                    if (i >= q.children[child_index].rb) {
                        ++child_index;
                    }
                }
            }
            for (/*nope*/; j<=q.rb; ++j) {
                pair_check(count_generated, pairs, i, j, SA, BWT, SID, DB, sentinal);
            }
        }
    }
    else {
        for (int i=q.lb; i<=q.rb; ++i) {
            for (int j=i+1; j<=q.rb; ++j) {
                pair_check(count_generated, pairs, i, j, SA, BWT, SID, DB, sentinal);
            }
        }
    }
}

inline static void print_array(
        const char * filename_,
        const int * const restrict array,
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len)
{
    int i;
    int j;
    FILE *f = NULL;
    const char *filename = filename_;

    f = fopen(filename, "w");
    if (NULL == f) {
        eprintf(stderr, "fopen(\"%s\") error: %s\n", filename, strerror(errno));
        exit(-1);
    }
    fprintf(f, " ");
    for (j=0; j<s2Len; ++j) {
        fprintf(f, "%4c", s2[j]);
    }
    fprintf(f, "\n");
    for (i=0; i<s1Len; ++i) {
        fprintf(f, "%c", s1[i]);
        for (j=0; j<s2Len; ++j) {
            fprintf(f, "%4d", array[i*s2Len + j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

inline static int self_score(
        const char * const restrict seq,
        int len,
        const parasail_matrix_t *matrix)
{
    int score = 0;
    for (int i=0; i<len; ++i) {
        unsigned char mapped = matrix->mapper[(unsigned char)seq[i]];
        score += matrix->matrix[matrix->size*mapped+mapped];
    }
    return score;
}

inline static void output_edges(
        FILE *fop,
        bool has_query,
        long sid_crossover,
        unsigned char *T,
        int AOL,
        int SIM,
        int OS,
        const parasail_matrix_t *matrix,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results,
        long long start,
        long long stop)
{
    unsigned long edge_count = 0;
    for (long long index=start; index<stop; ++index) {
        parasail_result_t *result = results[index];
        int i = vpairs[index].first;
        int j = vpairs[index].second;
        long i_beg = BEG[i];
        long i_end = END[i];
        long i_len = i_end-i_beg;
        long j_beg = BEG[j];
        long j_end = END[j];
        long j_len = j_end-j_beg;

        if (has_query) {
            i = i - sid_crossover;
        }

        if (parasail_result_is_saturated(result)) {
            if (has_query) {
                fprintf(stderr, "query %d and ref %d saturated\n", i, j);
            } else {
                fprintf(stderr, "seq %d and seq %d saturated\n", i, j);
            }
            continue;
        }

        int score = parasail_result_get_score(result);
        int matches = parasail_result_get_matches(result);
        int length = parasail_result_get_length(result);
        int self_score_ = 0;
        int max_len = 0;
        int i_self_score = self_score(
                (const char*)&T[i_beg], i_len, matrix);
        int j_self_score = self_score(
                (const char*)&T[j_beg], j_len, matrix);

        if (i_len > j_len) {
            max_len = i_len;
            self_score_ = i_self_score;
        }
        else {
            max_len = j_len;
            self_score_ = j_self_score;
        }

        if ((length * 100 >= AOL * int(max_len))
                && (matches * 100 >= SIM * length)
                && (score * 100 >= OS * self_score_)) {
            ++edge_count;
            fprintf(fop, "%d,%d,%f,%f,%f\n",
                    i, j,
                    1.0*length/max_len,
                    1.0*matches/length,
                    1.0*score/self_score_);
        }
    }

    fprintf(stdout, "%20s: %lu\n", "edges count", edge_count);
}

inline static void output_graph(
        FILE *fop,
        int which,
        unsigned char *T,
        int AOL,
        int SIM,
        int OS,
        const parasail_matrix_t *matrix,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results,
        vector<vector<pair<int,float> > > &graph,
        unsigned long &edge_count,
        long long start,
        long long stop)
{
    if (NULL == fop) {
        for (long long index=start; index<stop; ++index) {
            parasail_result_t *result = results[index];
            int i = vpairs[index].first;
            int j = vpairs[index].second;
            long i_beg = BEG[i];
            long i_end = END[i];
            long i_len = i_end-i_beg;
            long j_beg = BEG[j];
            long j_end = END[j];
            long j_len = j_end-j_beg;

            if (parasail_result_is_saturated(result)) {
                fprintf(stderr, "seq %d and seq %d saturated\n", i, j);
                continue;
            }

            int score = parasail_result_get_score(result);
            int matches = parasail_result_get_matches(result);
            int length = parasail_result_get_length(result);
            int self_score_ = 0;
            int max_len = 0;
            int i_self_score = self_score(
                    (const char*)&T[i_beg], i_len, matrix);
            int j_self_score = self_score(
                    (const char*)&T[j_beg], j_len, matrix);

            if (i_len > j_len) {
                max_len = i_len;
                self_score_ = i_self_score;
            }
            else {
                max_len = j_len;
                self_score_ = j_self_score;
            }

            if ((length * 100 >= AOL * int(max_len))
                    && (matches * 100 >= SIM * length)
                    && (score * 100 >= OS * self_score_)) {
                float value;
                ++edge_count;
                switch (which) {
                    case 0:
                        value = 1.0*length/max_len;
                        break;
                    case 1:
                        value = 1.0*matches/length;
                        break;
                    case 2:
                        value = 1.0*score/self_score_;
                        break;
                }
                graph[i].push_back(make_pair(j,value));
                graph[j].push_back(make_pair(i,value));
            }
        }
    }
    else {
        fprintf(fop, "%lu %lu 1\n", (unsigned long)graph.size(), edge_count);
        for (size_t i=0; i<graph.size(); ++i) {
            if (graph[i].size() > 0) {
                fprintf(fop, "%d %f", graph[i][0].first+1, graph[i][0].second);
                for (size_t j=1; j<graph[i].size(); ++j) {
                    fprintf(fop, "  %d %f", graph[i][j].first+1, graph[i][j].second);
                }
                fprintf(fop, "\n");
            } else {
                fprintf(fop, "\n");
            }
        }
        fprintf(stdout, "%20s: %lu\n", "edges count", edge_count);
    }
}

inline static void output_stats(
        FILE *fop,
        bool has_query,
        long sid_crossover,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results,
        long long start,
        long long stop)
{
    for (long long index=start; index<stop; ++index) {
        parasail_result_t *result = results[index];
        int i = vpairs[index].first;
        int j = vpairs[index].second;
        long i_beg = BEG[i];
        long i_end = END[i];
        long i_len = i_end-i_beg;
        long j_beg = BEG[j];
        long j_end = END[j];
        long j_len = j_end-j_beg;

        if (has_query) {
            i = i - sid_crossover;
        }

        if (parasail_result_is_saturated(result)) {
            if (has_query) {
                fprintf(stderr, "query %d and ref %d saturated\n", i, j);
            } else {
                fprintf(stderr, "seq %d and seq %d saturated\n", i, j);
            }
            continue;
        }

        eprintf(fop, "%d,%d,%ld,%ld,%d,%d,%d,%d,%d,%d\n",
                i,
                j,
                i_len,
                j_len,
                parasail_result_get_score(result),
                parasail_result_get_end_query(result),
                parasail_result_get_end_ref(result),
                parasail_result_get_matches(result),
                parasail_result_get_similar(result),
                parasail_result_get_length(result));
    }
}

inline static void output_basic(
        FILE *fop,
        bool has_query,
        long sid_crossover,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results,
        long long start,
        long long stop)
{
    for (long long index=start; index<stop; ++index) {
        parasail_result_t *result = results[index];
        int i = vpairs[index].first;
        int j = vpairs[index].second;
        long i_beg = BEG[i];
        long i_end = END[i];
        long i_len = i_end-i_beg;
        long j_beg = BEG[j];
        long j_end = END[j];
        long j_len = j_end-j_beg;

        if (has_query) {
            i = i - sid_crossover;
        }

        if (parasail_result_is_saturated(result)) {
            if (has_query) {
                fprintf(stderr, "query %d and ref %d saturated\n", i, j);
            } else {
                fprintf(stderr, "seq %d and seq %d saturated\n", i, j);
            }
            continue;
        }

        eprintf(fop, "%d,%d,%ld,%ld,%d,%d,%d\n",
                i,
                j,
                i_len,
                j_len,
                parasail_result_get_score(result),
                parasail_result_get_end_query(result),
                parasail_result_get_end_ref(result));
    }
}

inline static void output_emboss(
        FILE *fop,
        bool has_query,
        long sid_crossover,
        const parasail_matrix_t *matrix,
        const PairVec &vpairs,
        parasail_sequences_t *queries,
        parasail_sequences_t *sequences,
        const vector<parasail_result_t*> &results,
        long long start,
        long long stop)
{
    for (long long index=start; index<stop; ++index) {
        parasail_result_t *result = results[index];
        int i = vpairs[index].first;
        int j = vpairs[index].second;

        if (has_query) {
            i = i - sid_crossover;
            if (parasail_result_is_saturated(result)) {
                fprintf(stderr, "query %d (%s) and ref %d (%s) saturated\n",
                        i, queries->seqs[i].name.s,
                        j, sequences->seqs[j].name.s);
                continue;
            }
            parasail_traceback_generic_extra(
                    queries->seqs[i].seq.s,
                    queries->seqs[i].seq.l,
                    sequences->seqs[j].seq.s,
                    sequences->seqs[j].seq.l,
                    queries->seqs[i].name.s,
                    sequences->seqs[j].name.s,
                    matrix,
                    result,
                    '|', ':', '.',
                    50,
                    14,
                    1,
                    7,
                    fop);
        }
        else {
            if (parasail_result_is_saturated(result)) {
                fprintf(stderr, "seq %d (%s) and seq %d (%s) saturated\n",
                        i, sequences->seqs[i].name.s,
                        j, sequences->seqs[j].name.s);
                continue;
            }
            parasail_traceback_generic_extra(
                    sequences->seqs[i].seq.s,
                    sequences->seqs[i].seq.l,
                    sequences->seqs[j].seq.s,
                    sequences->seqs[j].seq.l,
                    sequences->seqs[i].name.s,
                    sequences->seqs[j].name.s,
                    matrix,
                    result,
                    '|', ':', '.',
                    50,
                    14,
                    1,
                    7,
                    fop);
        }
    }
}

inline static void output_ssw(
        FILE *fop,
        bool has_query,
        long sid_crossover,
        const parasail_matrix_t *matrix,
        const PairVec &vpairs,
        parasail_sequences_t *queries,
        parasail_sequences_t *sequences,
        const vector<parasail_result_t*> &results,
        long long start,
        long long stop)
{
    for (long long index=start; index<stop; ++index) {
        parasail_result_t *result = results[index];
        int i = vpairs[index].first;
        int j = vpairs[index].second;
        parasail_cigar_t *cigar = NULL;

        if (has_query) {
            i = i - sid_crossover;
            if (parasail_result_is_saturated(result)) {
                fprintf(stderr, "query %d (%s) and ref %d (%s) saturated\n",
                        i, queries->seqs[i].name.s,
                        j, sequences->seqs[j].name.s);
                continue;
            }
            fprintf(fop, "target_name: %s\n", sequences->seqs[j].name.s);
            fprintf(fop, "query_name: %s\n", queries->seqs[i].name.s);
            cigar = parasail_result_get_cigar(result,
                    queries->seqs[i].seq.s,
                    queries->seqs[i].seq.l,
                    sequences->seqs[j].seq.s,
                    sequences->seqs[j].seq.l,
                    matrix);
        }
        else {
            if (parasail_result_is_saturated(result)) {
                fprintf(stderr, "seq %d (%s) and seq %d (%s) saturated\n",
                        i, sequences->seqs[i].name.s,
                        j, sequences->seqs[j].name.s);
                continue;
            }
            fprintf(fop, "target_name: %s\n", sequences->seqs[j].name.s);
            fprintf(fop, "query_name: %s\n", sequences->seqs[i].name.s);
            cigar = parasail_result_get_cigar(result,
                    sequences->seqs[i].seq.s,
                    sequences->seqs[i].seq.l,
                    sequences->seqs[j].seq.s,
                    sequences->seqs[j].seq.l,
                    matrix);
        }

        fprintf(fop, "optimal_alignment_score: %d"
                "\tstrand: +"
                "\ttarget_begin: %d"
                "\ttarget_end: %d"
                "\tquery_begin: %d"
                "\tquery_end: %d\n",
                result->score,
                cigar->beg_ref+1,
                parasail_result_get_end_ref(result)+1,
                cigar->beg_query+1,
                parasail_result_get_end_query(result)+1);

        /* we only needed the cigar for beginning locations */
        parasail_cigar_free(cigar);

        if (has_query) {
            parasail_traceback_generic_extra(
                    queries->seqs[i].seq.s,
                    queries->seqs[i].seq.l,
                    sequences->seqs[j].seq.s,
                    sequences->seqs[j].seq.l,
                    "Query:",
                    "Target:",
                    matrix,
                    result,
                    '|', '*', '*',
                    60,
                    10,
                    0,
                    7,
                    fop);
        }
        else {
            parasail_traceback_generic_extra(
                    sequences->seqs[i].seq.s,
                    sequences->seqs[i].seq.l,
                    sequences->seqs[j].seq.s,
                    sequences->seqs[j].seq.l,
                    "Query:",
                    "Target:",
                    matrix,
                    result,
                    '|', '*', '*',
                    60,
                    10,
                    0,
                    7,
                    fop);
        }
    }
}

inline static void output_sam(
        FILE *fop,
        bool use_sam_header,
        bool has_query,
        long sid_crossover,
        const parasail_matrix_t *matrix,
        const PairVec &vpairs,
        parasail_sequences_t *queries,
        parasail_sequences_t *sequences,
        const vector<parasail_result_t*> &results,
        long long start,
        long long stop)
{
    if (use_sam_header && has_query && 0 == start) {
        fprintf(fop, "@HD\tVN:1.4\tSO:queryname\n");
        for (size_t index=0; index<sequences->l; ++index) {
            parasail_sequence_t ref_seq = sequences->seqs[index];
            fprintf(fop, "@SQ\tSN:%s\tLN:%d\n",
                    ref_seq.name.s, (int32_t)ref_seq.seq.l);
        }
    }
    for (long long index=start; index<stop; ++index) {
        parasail_result_t *result = results[index];
        int i = vpairs[index].first;
        int j = vpairs[index].second;
        parasail_sequence_t ref_seq;
        parasail_sequence_t read_seq;

        ref_seq = sequences->seqs[j];
        if (has_query) {
            i = i - sid_crossover;
            if (parasail_result_is_saturated(result)) {
                fprintf(stderr, "query %d (%s) and ref %d (%s) saturated\n",
                        i, queries->seqs[i].name.s,
                        j, sequences->seqs[j].name.s);
                continue;
            }
            read_seq = queries->seqs[i];
        }
        else {
            if (parasail_result_is_saturated(result)) {
                fprintf(stderr, "seq %d (%s) and seq %d (%s) saturated\n",
                        i, sequences->seqs[i].name.s,
                        j, sequences->seqs[j].name.s);
                continue;
            }
            read_seq = sequences->seqs[i];
        }

        fprintf(fop, "%s\t", read_seq.name.s);
        if (result->score == 0) {
            fprintf(fop, "4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
        }
        else {
            int32_t c = 0;
            int32_t length = 0;
            uint32_t mapq = 255; /* not available */
            parasail_cigar_t *cigar = NULL;
            uint32_t mismatch = 0;

            cigar = parasail_result_get_cigar(
                    result,
                    read_seq.seq.s, read_seq.seq.l,
                    ref_seq.seq.s, ref_seq.seq.l,
                    matrix);

            fprintf(fop, "0\t");
            fprintf(fop, "%s\t%d\t%d\t",
                    ref_seq.name.s, cigar->beg_ref + 1, mapq);
            if (parasail_result_is_sw(result)) {
                if (cigar->beg_query > 0) {
                    fprintf(fop, "%dS", cigar->beg_query);
                }
            }
            for (c=0; c<cigar->len; ++c) {
                char letter = parasail_cigar_decode_op(cigar->seq[c]);
                uint32_t length = parasail_cigar_decode_len(cigar->seq[c]);
                fprintf(fop, "%lu%c", (unsigned long)length, letter);
                if ('X' == letter || 'I' == letter || 'D' == letter) {
                    mismatch += length;
                }
            }

            length = read_seq.seq.l - result->end_query - 1;
            if (parasail_result_is_sw(result)) {
                if (length > 0) {
                    fprintf(fop, "%dS", length);
                }
            }
            fprintf(fop, "\t*\t0\t0\t");
            fprintf(fop, "%s", read_seq.seq.s);
            fprintf(fop, "\t");
            if (read_seq.qual.s) {
                fprintf (fop, "%s", read_seq.qual.s);
            }
            else {
                fprintf(fop, "*");
            }
            fprintf(fop, "\tAS:i:%d", result->score);
            fprintf(fop,"\tNM:i:%d\t", mismatch);
            fprintf(fop, "\n");

            parasail_cigar_free(cigar);
        }
    }
}

inline static void output_trace(
        FILE *fop,
        bool has_query,
        long sid_crossover,
        unsigned char *T,
        const parasail_matrix_t *matrix,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results,
        long long start,
        long long stop)
{
    for (long long index=start; index<stop; ++index) {
        parasail_result_t *result = results[index];
        int i = vpairs[index].first;
        int j = vpairs[index].second;
        long i_beg = BEG[i];
        long i_end = END[i];
        long i_len = i_end-i_beg;
        long j_beg = BEG[j];
        long j_end = END[j];
        long j_len = j_end-j_beg;
        parasail_cigar_t *cigar = NULL;
        char *cigar_string = NULL;

        if (has_query) {
            i = i - sid_crossover;
        }

        if (parasail_result_is_saturated(result)) {
            if (has_query) {
                fprintf(stderr, "query %d and ref %d saturated\n", i, j);
            } else {
                fprintf(stderr, "seq %d and seq %d saturated\n", i, j);
            }
            continue;
        }

        cigar = parasail_result_get_cigar(
                result,
                (const char*)&T[i_beg], i_len,
                (const char*)&T[j_beg], j_len,
                matrix);
        cigar_string = parasail_cigar_decode(cigar);
        eprintf(fop, "%d,%d,%ld,%ld,%d,%d,%d,%s\n",
                i,
                j,
                i_len,
                j_len,
                parasail_result_get_score(result),
                parasail_result_get_end_query(result),
                parasail_result_get_end_ref(result),
                cigar_string);
        parasail_cigar_free(cigar);
        free(cigar_string);
    }
}

inline static void output_tables(
        bool has_query,
        long sid_crossover,
        unsigned char *T,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results,
        long long start,
        long long stop)
{
    for (long long index=start; index<stop; ++index) {
        parasail_result_t *result = results[index];
        int i = vpairs[index].first;
        int j = vpairs[index].second;
        long i_beg = BEG[i];
        long i_end = END[i];
        long i_len = i_end-i_beg;
        long j_beg = BEG[j];
        long j_end = END[j];
        long j_len = j_end-j_beg;
        int *table = parasail_result_get_score_table(result);

        if (has_query) {
            i = i - sid_crossover;
        }

        char filename[256] = {'\0'};
        sprintf(filename, "parasail_%d_%d.txt", i, j);
        print_array(filename, table,
                (const char*)&T[i_beg], i_len,
                (const char*)&T[j_beg], j_len);
    }
}

inline static void output(
        bool is_stats,
        bool is_table,
        bool is_trace,
        bool edge_output,
        bool use_emboss_format,
        bool use_ssw_format,
        bool use_sam_format,
        bool use_sam_header,
        FILE *fop,
        bool has_query,
        long sid_crossover,
        unsigned char *T,
        int AOL,
        int SIM,
        int OS,
        const parasail_matrix_t *matrix,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        parasail_sequences_t *queries,
        parasail_sequences_t *sequences,
        const vector<parasail_result_t*> &results,
        long long start,
        long long stop)
{
    if (is_stats) {
        if (edge_output) {
            output_edges(fop, has_query, sid_crossover, T, AOL, SIM, OS, matrix, BEG, END, vpairs, results, start, stop);
        }
        else {
            output_stats(fop, has_query, sid_crossover, BEG, END, vpairs, results, start, stop);
        }
    }
    else if (is_trace) {
        if (use_emboss_format) {
            output_emboss(fop, has_query, sid_crossover, matrix, vpairs, queries, sequences, results, start, stop);
        }
        else if (use_ssw_format) {
            output_ssw(fop, has_query, sid_crossover, matrix, vpairs, queries, sequences, results, start, stop);
        }
        else if (use_sam_format) {
            output_sam(fop, use_sam_header, has_query, sid_crossover, matrix, vpairs, queries, sequences, results, start, stop);
        }
        else {
            output_trace(fop, has_query, sid_crossover, T, matrix, BEG, END, vpairs, results, start, stop);
        }
    }
    else {
        output_basic(fop, has_query, sid_crossover, BEG, END, vpairs, results, start, stop);
    }
    if (is_table) {
        output_tables(has_query, sid_crossover, T, BEG, END, vpairs, results, start, stop);
    }
}

inline static size_t parse_bytes(const char *value)
{
    size_t multiplier = 0;
    double base = 0;
    string unit;
    istringstream iss(value);

    iss >> base;
    if (!iss) {
        eprintf(stderr, "could not parse memory value\n");
        return 0;
    }
    iss >> unit;
    if (!iss || unit.empty()) {
        eprintf(stderr, "could not parse memory unit, default to bytes\n");
        unit = "b";
    }
    transform(unit.begin(), unit.end(), unit.begin(), ::tolower);

    if ("b" == unit
            || "byte" == unit
            || "bytes" == unit) {
        multiplier = 1ULL;
    }
    else if ("kb" == unit
            || "kbyte" == unit
            || "kbytes" == unit
            || "kilobyte" == unit
            || "kilobytes" == unit) {
        multiplier = 1000ULL;
    }
    else if ("mb" == unit
            || "mbyte" == unit
            || "mbytes" == unit
            || "megabyte" == unit
            || "megabytes" == unit) {
        multiplier = 1000ULL * 1000ULL;
    }
    else if ("gb" == unit
            || "gbyte" == unit
            || "gbytes" == unit
            || "gigabyte" == unit
            || "gigabytes" == unit) {
        multiplier = 1000ULL * 1000ULL * 1000ULL;
    }
    else {
        eprintf(stderr, "could not parse byte unit\n");
        return 0;
    }

    return static_cast<size_t>(multiplier*base);
}

#ifdef HAVE_SETUNHANDLEDEXCEPTIONFILTER
LONG WINAPI windows_exception_handler(EXCEPTION_POINTERS * ExceptionInfo)
{
    switch (ExceptionInfo->ExceptionRecord->ExceptionCode)
    {
    case EXCEPTION_ACCESS_VIOLATION:
        fputs("Error: EXCEPTION_ACCESS_VIOLATION\n", stderr);
        break;
    case EXCEPTION_ARRAY_BOUNDS_EXCEEDED:
        fputs("Error: EXCEPTION_ARRAY_BOUNDS_EXCEEDED\n", stderr);
        break;
    case EXCEPTION_BREAKPOINT:
        fputs("Error: EXCEPTION_BREAKPOINT\n", stderr);
        break;
    case EXCEPTION_DATATYPE_MISALIGNMENT:
        fputs("Error: EXCEPTION_DATATYPE_MISALIGNMENT\n", stderr);
        break;
    case EXCEPTION_FLT_DENORMAL_OPERAND:
        fputs("Error: EXCEPTION_FLT_DENORMAL_OPERAND\n", stderr);
        break;
    case EXCEPTION_FLT_DIVIDE_BY_ZERO:
        fputs("Error: EXCEPTION_FLT_DIVIDE_BY_ZERO\n", stderr);
        break;
    case EXCEPTION_FLT_INEXACT_RESULT:
        fputs("Error: EXCEPTION_FLT_INEXACT_RESULT\n", stderr);
        break;
    case EXCEPTION_FLT_INVALID_OPERATION:
        fputs("Error: EXCEPTION_FLT_INVALID_OPERATION\n", stderr);
        break;
    case EXCEPTION_FLT_OVERFLOW:
        fputs("Error: EXCEPTION_FLT_OVERFLOW\n", stderr);
        break;
    case EXCEPTION_FLT_STACK_CHECK:
        fputs("Error: EXCEPTION_FLT_STACK_CHECK\n", stderr);
        break;
    case EXCEPTION_FLT_UNDERFLOW:
        fputs("Error: EXCEPTION_FLT_UNDERFLOW\n", stderr);
        break;
    case EXCEPTION_ILLEGAL_INSTRUCTION:
        fputs("Error: EXCEPTION_ILLEGAL_INSTRUCTION\n", stderr);
        break;
    case EXCEPTION_IN_PAGE_ERROR:
        fputs("Error: EXCEPTION_IN_PAGE_ERROR\n", stderr);
        break;
    case EXCEPTION_INT_DIVIDE_BY_ZERO:
        fputs("Error: EXCEPTION_INT_DIVIDE_BY_ZERO\n", stderr);
        break;
    case EXCEPTION_INT_OVERFLOW:
        fputs("Error: EXCEPTION_INT_OVERFLOW\n", stderr);
        break;
    case EXCEPTION_INVALID_DISPOSITION:
        fputs("Error: EXCEPTION_INVALID_DISPOSITION\n", stderr);
        break;
    case EXCEPTION_NONCONTINUABLE_EXCEPTION:
        fputs("Error: EXCEPTION_NONCONTINUABLE_EXCEPTION\n", stderr);
        break;
    case EXCEPTION_PRIV_INSTRUCTION:
        fputs("Error: EXCEPTION_PRIV_INSTRUCTION\n", stderr);
        break;
    case EXCEPTION_SINGLE_STEP:
        fputs("Error: EXCEPTION_SINGLE_STEP\n", stderr);
        break;
    case EXCEPTION_STACK_OVERFLOW:
        fputs("Error: EXCEPTION_STACK_OVERFLOW\n", stderr);
        break;
    default:
        fputs("Error: Unrecognized Exception\n", stderr);
        break;
    }
    fflush(stderr);
    return EXCEPTION_EXECUTE_HANDLER;
}
void set_signal_handler()
{
    SetUnhandledExceptionFilter(windows_exception_handler);
}
#else
void posix_signal_handler(int sig, siginfo_t *siginfo, void *context)
{
    (void)context;
    switch (sig)
    {
    case SIGSEGV:
        fputs("Caught SIGSEGV: Segmentation Fault\n", stderr);
        break;
    case SIGINT:
        fputs("Caught SIGINT: Interactive attention signal, (usually ctrl+c)\n", stderr);
        break;
    case SIGFPE:
        switch (siginfo->si_code)
        {
        case FPE_INTDIV:
            fputs("Caught SIGFPE: (integer divide by zero)\n", stderr);
            break;
        case FPE_INTOVF:
            fputs("Caught SIGFPE: (integer overflow)\n", stderr);
            break;
        case FPE_FLTDIV:
            fputs("Caught SIGFPE: (floating-point divide by zero)\n", stderr);
            break;
        case FPE_FLTOVF:
            fputs("Caught SIGFPE: (floating-point overflow)\n", stderr);
            break;
        case FPE_FLTUND:
            fputs("Caught SIGFPE: (floating-point underflow)\n", stderr);
            break;
        case FPE_FLTRES:
            fputs("Caught SIGFPE: (floating-point inexact result)\n", stderr);
            break;
        case FPE_FLTINV:
            fputs("Caught SIGFPE: (floating-point invalid operation)\n", stderr);
            break;
        case FPE_FLTSUB:
            fputs("Caught SIGFPE: (subscript out of range)\n", stderr);
            break;
        default:
            fputs("Caught SIGFPE: Arithmetic Exception\n", stderr);
            break;
        }
    case SIGILL:
        switch (siginfo->si_code)
        {
        case ILL_ILLOPC:
            fputs("Caught SIGILL: (illegal opcode)\n", stderr);
            break;
        case ILL_ILLOPN:
            fputs("Caught SIGILL: (illegal operand)\n", stderr);
            break;
        case ILL_ILLADR:
            fputs("Caught SIGILL: (illegal addressing mode)\n", stderr);
            break;
        case ILL_ILLTRP:
            fputs("Caught SIGILL: (illegal trap)\n", stderr);
            break;
        case ILL_PRVOPC:
            fputs("Caught SIGILL: (privileged opcode)\n", stderr);
            break;
        case ILL_PRVREG:
            fputs("Caught SIGILL: (privileged register)\n", stderr);
            break;
        case ILL_COPROC:
            fputs("Caught SIGILL: (coprocessor error)\n", stderr);
            break;
        case ILL_BADSTK:
            fputs("Caught SIGILL: (internal stack error)\n", stderr);
            break;
        default:
            fputs("Caught SIGILL: Illegal Instruction\n", stderr);
            break;
        }
        break;
    case SIGTERM:
        fputs("Caught SIGTERM: a termination request was sent to the program\n", stderr);
        break;
    case SIGABRT:
        fputs("Caught SIGABRT: usually caused by an abort() or assert()\n", stderr);
        break;
    default:
        break;
    }
    fflush(stderr);
    _exit(EXIT_FAILURE);
}
void set_signal_handler()
{
    /* register our signal handlers */
    {
        struct sigaction sig_action = {};
        sig_action.sa_sigaction = posix_signal_handler;
        sigemptyset(&sig_action.sa_mask);
        sig_action.sa_flags = SA_SIGINFO;
        if (sigaction(SIGSEGV, &sig_action, NULL) != 0) { err(1, "sigaction"); }
        if (sigaction(SIGFPE, &sig_action, NULL) != 0) { err(1, "sigaction"); }
        if (sigaction(SIGINT, &sig_action, NULL) != 0) { err(1, "sigaction"); }
        if (sigaction(SIGILL, &sig_action, NULL) != 0) { err(1, "sigaction"); }
        if (sigaction(SIGTERM, &sig_action, NULL) != 0) { err(1, "sigaction"); }
        if (sigaction(SIGABRT, &sig_action, NULL) != 0) { err(1, "sigaction"); }
    }
}
#endif

