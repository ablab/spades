#include "config.h"

#include <ctype.h>
#include <stdlib.h>

#if defined(_MSC_VER)
#include "wingetopt/src/getopt.h"
#else
#include <unistd.h>
#endif

#include "parasail/parasail.h"
#include "parasail/parasail/io.h"
#include "ssw.h"


/* This table is used to transform amino acid letters into numbers. */
int8_t aa_table[128] = {
	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
	14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
	23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
	14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23
};

/* This table is used to transform nucleotide letters into numbers. */
int8_t nt_table[128] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


static void print_help(const char *progname, int status) {
    fprintf(stderr, "\nusage: %s "
            "[-e gap_extend] "
            "[-o gap_open] "
            "[-m matrix] "
            "[-d] "
            "[-M match] "
            "[-X mismatch] "
            "-f file "
            "\n\n",
            progname);
    fprintf(stderr, "Defaults:\n"
            " gap_extend: 1, must be >= 0\n"
            "   gap_open: 10, must be >= 0\n"
            "     matrix: blosum62\n"
            "         -d: if present, assume DNA alphabet\n"
            "      match: 1, must be >= 0\n"
            "   mismatch: 0, must be >= 0\n"
            "       file: no default, must be in FASTA format\n"
            );
    exit(status);
}

static void print_cigar(const uint32_t *cigar, int32_t cigarLen)
{
    int32_t i;

    for (i=0; i<cigarLen; ++i) {
        char c = parasail_cigar_decode_op(cigar[i]);
        uint32_t l = parasail_cigar_decode_len(cigar[i]);
        printf("%u%c", l, c);
    }
    printf("\n");
}


int main(int argc, char **argv)
{
    const char *fname = NULL;
    parasail_sequences_t *sequences = NULL;
    double start = 0;
    double finish = 0;
    int i = 0;
    size_t j = 0;
    int c = 0;
    const char *matrixname = NULL;
    const parasail_matrix_t *matrix = NULL;
    int8_t ssw_matrix[24*24];
    int gap_open = 10;
    int gap_extend = 1;
    int match = 1;
    int mismatch = 0;
    int use_dna = 0;
    int8_t *table = aa_table;
    int32_t n = 24;
    const char *progname = "ssw_vs_pssw";
    int do_time = 0;

    /* Check arguments. */
    while ((c = getopt(argc, argv, "de:f:hm:M:o:tX:")) != -1) {
        switch (c) {
            case 'd':
                use_dna = 1;
                break;
            case 'e':
                gap_extend = atoi(optarg);
                if (gap_extend < 0) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 'f':
                fname = optarg;
                break;
            case 'h':
                print_help(progname, EXIT_FAILURE);
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
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 't':
                do_time = 1;
                break;
            case 'X':
                mismatch = atoi(optarg);
                if (mismatch < 0) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case '?':
                if ( optopt == 'e'
                        || optopt == 'm'
                        || optopt == 'M'
                        || optopt == 'o'
                        || optopt == 'X'
                   ) {
                    fprintf(stderr,
                            "Option -%c requires an argument.\n",
                            optopt);
                }
                else if (isprint(optopt)) {
                    fprintf(stderr, "Unknown option `-%c'.\n",
                            optopt);
                }
                else {
                    fprintf(stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                }
                print_help(progname, EXIT_FAILURE);
                break;
            default:
                fprintf(stderr, "default case in getopt\n");
                exit(EXIT_FAILURE);
        }
    }

    /* select the substitution matrix */
    if (NULL != matrixname && use_dna) {
        fprintf(stderr, "Cannot specify matrix name for DNA alignments.\n");
        exit(EXIT_FAILURE);
    }
    if (use_dna) {
        matrix = parasail_matrix_create("ACGT", match, -mismatch);
		table = nt_table;
        n = 5;
    }
    else {
        if (NULL == matrixname) {
            matrixname = "blosum62";
        }
        matrix = parasail_matrix_lookup(matrixname);
        if (NULL == matrix) {
            fprintf(stderr, "Specified substitution matrix not found.\n");
            exit(EXIT_FAILURE);
        }
    }
    /* create the ssw matrix */
    for (i=0; i<matrix->size*matrix->size; ++i) {
        ssw_matrix[i] = matrix->matrix[i];
    }

    if (fname == NULL) {
        fprintf(stderr, "missing input file\n");
        print_help(progname, EXIT_FAILURE);
    }

    /* print the parameters for reference */
    fprintf(stdout,
            "%20s: %d\n"
            "%20s: %d\n"
            "%20s: %s\n"
            "%20s: %s\n",
            "gap_extend", gap_extend,
            "gap_open", gap_open,
            "matrix", matrixname,
            "file", fname);

    sequences = parasail_sequences_from_file(fname);

    if (do_time) {
        /* ssw timing */
        start = parasail_time();
        for (j=0; j<sequences->l; ++j) {
            size_t k = 0;
            char *readSeq = NULL;
            int32_t readLen = 0;
            int32_t maskLen = 0;
            s_profile *p = NULL;
            int8_t *num = NULL;
            int32_t m = 0;

            /* prep ssw */
            readSeq = sequences->seqs[j].seq.s;
            readLen = sequences->seqs[j].seq.l;
            maskLen = readLen / 2;
            num = (int8_t*)malloc(readLen);
            for (m=0; m<readLen; ++m) num[m] = table[(int)readSeq[m]];
            p = ssw_init(num, readLen, ssw_matrix, n, 2);

            for (k=j+1; k<sequences->l; ++k) {
                int8_t *ref_num = NULL;
                s_align *result = NULL;
                char *refSeq = NULL;
                int32_t refLen = 0;
                int32_t filter = 0;
                int8_t flag = 2;

                refSeq = sequences->seqs[k].seq.s;
                refLen = sequences->seqs[k].seq.l;
                ref_num = (int8_t*)malloc(refLen);
                for (m=0; m<refLen; ++m) ref_num[m] = table[(int)refSeq[m]];
                result = ssw_align(p, ref_num, refLen, gap_open, gap_extend, flag, filter, 0, maskLen);
                align_destroy(result);
                free(ref_num);
            }

            init_destroy(p);
            free(num);
        }
        finish = parasail_time();
        printf("  ssw: %f\n", finish-start);
        /* pssw timing */
        start = parasail_time();
        for (j=0; j<sequences->l; ++j) {
            size_t k = 0;
            char *readSeq = NULL;
            int32_t readLen = 0;

            /* prep ssw */
            readSeq = sequences->seqs[j].seq.s;
            readLen = sequences->seqs[j].seq.l;

            for (k=j+1; k<sequences->l; ++k) {
                char *refSeq = NULL;
                int32_t refLen = 0;
                parasail_result_ssw_t *result = NULL;

                refSeq = sequences->seqs[k].seq.s;
                refLen = sequences->seqs[k].seq.l;
                /* now pssw */
                result = parasail_ssw(readSeq, readLen, refSeq, refLen, gap_open, gap_extend, matrix);
                parasail_result_ssw_free(result);
            }
        }
        finish = parasail_time();
        printf(" pssw: %f\n", finish-start);
        /* trace timing */
        start = parasail_time();
        for (j=0; j<sequences->l; ++j) {
            size_t k = 0;
            char *readSeq = NULL;
            int32_t readLen = 0;

            readSeq = sequences->seqs[j].seq.s;
            readLen = sequences->seqs[j].seq.l;

            for (k=j+1; k<sequences->l; ++k) {
                char *refSeq = NULL;
                int32_t refLen = 0;
                parasail_result_t *result = NULL;
                parasail_cigar_t *cigar = NULL;

                refSeq = sequences->seqs[k].seq.s;
                refLen = sequences->seqs[k].seq.l;
                result = parasail_sw_trace_striped_8(readSeq, readLen, refSeq, refLen, gap_open, gap_extend, matrix);
                if (parasail_result_is_saturated(result)) {
                    parasail_result_free(result);
                    result = parasail_sw_trace_striped_16(readSeq, readLen, refSeq, refLen, gap_open, gap_extend, matrix);
                }
                cigar = parasail_result_get_cigar(result,
                        readSeq, readLen, refSeq, refLen, matrix);
                parasail_cigar_free(cigar);
                parasail_result_free(result);
            }
        }
        finish = parasail_time();
        printf("trace: %f\n", finish-start);
    }
    else {
        for (j=0; j<sequences->l; ++j) {
            size_t k = 0;
            char *readSeq = NULL;
            int32_t readLen = 0;
            int32_t maskLen = 0;
            s_profile *p = NULL;
            int8_t *num = NULL;
            int32_t m = 0;

            /* prep ssw */
            readSeq = sequences->seqs[j].seq.s;
            readLen = sequences->seqs[j].seq.l;
            maskLen = readLen / 2;
            num = (int8_t*)malloc(readLen);
            for (m=0; m<readLen; ++m) num[m] = table[(int)readSeq[m]];
            p = ssw_init(num, readLen, ssw_matrix, n, 2);

            for (k=j+1; k<sequences->l; ++k) {
                int8_t *ref_num = NULL;
                s_align *result = NULL;
                char *refSeq = NULL;
                int32_t refLen = 0;
                int32_t filter = 0;
                int8_t flag = 2;
                parasail_result_ssw_t *presult = NULL;
                parasail_result_t *tresult = NULL;
                parasail_cigar_t *cigar = NULL;

                printf("read %lu ref %lu\n", (long unsigned)j, (long unsigned)k);

                refSeq = sequences->seqs[k].seq.s;
                refLen = sequences->seqs[k].seq.l;
                ref_num = (int8_t*)malloc(refLen);
                for (m=0; m<refLen; ++m) ref_num[m] = table[(int)refSeq[m]];
                result = ssw_align(p, ref_num, refLen, gap_open, gap_extend, flag, filter, 0, maskLen);

                printf("--SSW--\n");
                printf("      score: %d\n", result->score1);
                printf(" ref_begin1: %d\n", result->ref_begin1);
                printf("   ref_end1: %d\n", result->ref_end1);
                printf("read_begin1: %d\n", result->read_begin1);
                printf("  read_end1: %d\n", result->read_end1);
                printf("   cigarLen: %d\n", result->cigarLen);
                printf("      cigar: ");
                print_cigar(result->cigar, result->cigarLen);

                align_destroy(result);
                free(ref_num);

                /* now pssw */
                presult = parasail_ssw(readSeq, readLen, refSeq, refLen, gap_open, gap_extend, matrix);

                printf("--PSSW--\n");
                printf("      score: %d\n", presult->score1);
                printf(" ref_begin1: %d\n", presult->ref_begin1);
                printf("   ref_end1: %d\n", presult->ref_end1);
                printf("read_begin1: %d\n", presult->read_begin1);
                printf("  read_end1: %d\n", presult->read_end1);
                printf("   cigarLen: %d\n", presult->cigarLen);
                printf("      cigar: ");
                print_cigar(presult->cigar, presult->cigarLen);

                parasail_result_ssw_free(presult);

                /* now trace */
                tresult = parasail_sw_trace_striped_sat(readSeq, readLen, refSeq, refLen, gap_open, gap_extend, matrix);
                cigar = parasail_result_get_cigar(tresult,
                        readSeq, readLen, refSeq, refLen, matrix);

                printf("--TRACE--\n");
                printf("      score: %d\n", tresult->score);
                printf("   ref_end1: %d\n", tresult->end_ref);
                printf("  read_end1: %d\n", tresult->end_query);
                printf("   cigarLen: %d\n", cigar->len);
                printf("      cigar: ");
                print_cigar(cigar->seq, cigar->len);

                parasail_cigar_free(cigar);
                parasail_result_free(tresult);
            }

            init_destroy(p);
            free(num);
        }
    }

    parasail_sequences_free(sequences);

    return 0;
}
