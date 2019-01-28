#include "config.h"

#include <assert.h>
#include <ctype.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#if defined(_MSC_VER)
#include "wingetopt/src/getopt.h"
#else
#include <unistd.h>
#endif

#include "parasail/parasail.h"
#include "parasail/parasail/io.h"
#include "parasail/parasail/memory.h"

static char* strdup_rpl(const char *str)
{
    size_t l = strlen(str);
    char *s = (char*)calloc(l + 1, sizeof(char));
    if(!s) {
        return NULL;
    }
    (void)memcpy(s, str, l);
    return s;
}

int main(int argc, char **argv)
{
    long seqA_index = 0;
    long seqB_index = 1;
    const char *seqA = NULL;
    const char *seqB = NULL;
    int lena = 0;
    int lenb = 0;
    parasail_result_t *result = NULL;
    int c = 0;
    char *filename = NULL;
    parasail_sequences_t *sequences = NULL;
    unsigned long seq_count = 0;
    char *endptr = NULL;
    const char *funcname = "sw_trace";
    parasail_function_t *function = NULL;
    parasail_pfunction_t *pfunction = NULL;
    parasail_pcreator_t *pcreator = NULL;
    parasail_profile_t *profile = NULL;
    char *matrixname = NULL;
    const parasail_matrix_t *matrix = NULL;
    int open = 10;
    int extend = 1;
    int match = 1;
    int mismatch = 0;
    int use_dna = 0;
    parasail_cigar_t *cigar = NULL;
    parasail_cigar_t *cigar2 = NULL;
    char *cigar_string = NULL;
    parasail_traceback_t *traceback = NULL;

    while ((c = getopt(argc, argv, "a:x:y:f:m:o:e:dX:Y:")) != -1) {
        switch (c) {
            case 'a':
                funcname = optarg;
                break;
            case 'x':
                errno = 0;
                seqA_index = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol seqA_index");
                    fprintf(stderr, "invalid seqA index\n");
                    exit(1);
                }
                break;
            case 'y':
                errno = 0;
                seqB_index = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol seqB_index");
                    fprintf(stderr, "invalid seqB index\n");
                    exit(1);
                }
                break;
            case 'f':
                filename = optarg;
                break;
            case 'd':
                use_dna = 1;
                break;
            case 'X':
                match = atoi(optarg);
                if (match < 0) {
                    fprintf(stderr, "match must be >= 0\n");
                    exit(1);
                }
                break;
            case 'Y':
                mismatch = atoi(optarg);
                if (mismatch < 0) {
                    fprintf(stderr, "mismatch must be >= 0\n");
                    exit(1);
                }
                break;
            case 'm':
                matrixname = optarg;
                break;
            case 'o':
                errno = 0;
                open = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol open");
                    exit(1);
                }
                break;
            case 'e':
                errno = 0;
                extend = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol extend");
                    exit(1);
                }
                break;
            case '?':
                if (optopt == 'a'
                        || optopt == 'b'
                        || optopt == 'f'
                        || optopt == 'n') {
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
                exit(1);
            default:
                fprintf(stderr, "default case in getopt\n");
                exit(1);
        }
    }

    if (filename) {
        sequences = parasail_sequences_from_file(filename);
        seq_count = sequences->l;
    }
    else {
        fprintf(stderr, "no filename specified\n");
        exit(1);
    }

    /* select the substitution matrix */
    if (NULL != matrixname && use_dna) {
        fprintf(stderr, "Cannot specify matrix name for DNA alignments.\n");
        exit(EXIT_FAILURE);
    }
    if (use_dna) {
        matrix = parasail_matrix_create("ACGT", match, -mismatch);
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


    if ((unsigned long)seqA_index >= seq_count) {
        fprintf(stderr, "seqA index out of bounds\n");
        exit(1);
    }

    if ((unsigned long)seqB_index >= seq_count) {
        fprintf(stderr, "seqB index out of bounds\n");
        exit(1);
    }

    /* make sure we were given a function that traces */
    if (NULL ==  strstr(funcname, "trace")) {
        fprintf(stderr, "function must enable backtraces\n");
        exit(1);
    }

    /* select the function */
    if (funcname) {
        if (NULL != strstr(funcname, "profile")) {
            pfunction = parasail_lookup_pfunction(funcname);
            if (NULL == pfunction) {
                fprintf(stderr, "Specified profile function not found.\n");
                exit(EXIT_FAILURE);
            }
            pcreator = parasail_lookup_pcreator(funcname);
            if (NULL == pcreator) {
                fprintf(stderr, "Specified profile creator function not found.\n");
                exit(EXIT_FAILURE);
            }

        }
        else {
            function = parasail_lookup_function(funcname);
            if (NULL == function) {
                fprintf(stderr, "Specified function not found.\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    else {
        fprintf(stderr, "No alignment function specified.\n");
        exit(EXIT_FAILURE);
    }

    seqA = sequences->seqs[seqA_index].seq.s;
    seqB = sequences->seqs[seqB_index].seq.s;
    lena = sequences->seqs[seqA_index].seq.l;
    lenb = sequences->seqs[seqB_index].seq.l;

    printf("        file: %s\n", filename);
    printf("    function: %s\n", funcname);
    printf("      matrix: %s\n", matrixname);
    printf("    gap open: %d\n", open);
    printf("  gap extend: %d\n", extend);
    printf("    seq pair: %lu,%lu\n", seqA_index, seqB_index);
    printf("query length: %d\n", lena);
    printf("   db length: %d\n", lenb);

    if (pfunction) {
        profile = pcreator(seqA, lena, matrix);
        result = pfunction(profile, seqB, lenb, open, extend);
        parasail_profile_free(profile);
    }
    else {
        result = function(seqA, lena, seqB, lenb, open, extend, matrix);
    }
    assert(parasail_result_is_trace(result));
    assert(result->trace->trace_table);

    /* do the traceback */
    parasail_traceback_generic(seqA, lena, seqB, lenb,
            sequences->seqs[seqA_index].name.s,
            sequences->seqs[seqB_index].name.s,
            matrix, result,
            '|', ':', '.', 50, 17, 1);

    parasail_traceback_generic(seqA, lena, seqB, lenb,
            "Query:",
            "Target:",
            matrix, result,
            '|', '*', '*', 60, 7, 0);

    printf("Score:          %d\n", result->score);
    printf("end_query:      %d\n", result->end_query);
    printf("end_ref:        %d\n", result->end_ref);

    /* test new traceback string functions */
    traceback = parasail_result_get_traceback(result, seqA, lena, seqB, lenb, matrix, '|', ':', '.');
    printf("\nTraceback string function\n");
    printf("%s\n", traceback->query);
    printf("%s\n", traceback->comp);
    printf("%s\n", traceback->ref);
    printf("\n");
    parasail_traceback_free(traceback);

    /* do the cigar */
    cigar = parasail_result_get_cigar(result, seqA, lena, seqB, lenb, matrix);
    printf("cigar uint32_t: '");
    for (c=0; c<cigar->len; ++c) {
        printf("%u", cigar->seq[c]);
    }
    printf("'\n");
    printf("cigar string:   '");
    for (c=0; c<cigar->len; ++c) {
        uint32_t len = parasail_cigar_decode_len(cigar->seq[c]);
        char op =  parasail_cigar_decode_op(cigar->seq[c]);
        printf("%u%c", len, op);
    }
    printf("'\n");

    printf("beg_query:      %d\n", cigar->beg_query);
    printf("beg_ref:        %d\n", cigar->beg_ref);

    /* translate cigar into traceback */
    {
        uint32_t i;
        uint32_t j;
        char *_a = malloc(lena+lenb+1);
        char *_b = malloc(lena+lenb+1);
        char *_x = malloc(lena+lenb+1);
        char *a = _a;
        char *b = _b;
        char *x = _x;
        size_t ao = cigar->beg_query;
        size_t bo = cigar->beg_ref;
        uint32_t size = 0;
        for (c=0; c<cigar->len; ++c) {
            uint32_t len = parasail_cigar_decode_len(cigar->seq[c]);
            char op =  parasail_cigar_decode_op(cigar->seq[c]);
            size += len;
            for (i=0; i<len; ++i) {
                if ('=' == op) {
                    *(a++) = sequences->seqs[seqA_index].seq.s[ao++];
                    *(b++) = sequences->seqs[seqB_index].seq.s[bo++];
                    *(x++) = '|';
                }
                else if ('X' == op) {
                    *(a++) = sequences->seqs[seqA_index].seq.s[ao++];
                    *(b++) = sequences->seqs[seqB_index].seq.s[bo++];
                    *(x++) = '*';
                }
                else if ('I' == op) {
                    *(a++) = sequences->seqs[seqA_index].seq.s[ao++];
                    *(b++) = '-';
                    *(x++) = ' ';
                }
                else if ('D' == op) {
                    *(a++) = '-';
                    *(b++) = sequences->seqs[seqB_index].seq.s[bo++];
                    *(x++) = ' ';
                }
            }
        }
        *a = '\0';
        *b = '\0';
        *x = '\0';
        a = _a;
        b = _b;
        x = _x;
        for (i=0; i<size; i+=60) {
            printf("\n");
            for (j=0; j<size&&j<60&&*a!='\0'; ++j) {
                printf("%c", *(a++));
            }
            printf("\n");
            for (j=0; j<size&&j<60&&*x!='\0'; ++j) {
                printf("%c", *(x++));
            }
            printf("\n");
            for (j=0; j<size&&j<60&&*b!='\0'; ++j) {
                printf("%c", *(b++));
            }
            printf("\n");
        }
        printf("\n");
        free(_a);
        free(_b);
        free(_x);
    }

    /* test the encode/decode functionality */
    cigar_string = parasail_cigar_decode(cigar);
    printf("cigar string:   '%s'\n", cigar_string);
    cigar2 = parasail_cigar_encode_string(cigar_string);
    printf("cigar2 uint32_t:'");
    for (c=0; c<cigar2->len; ++c) {
        printf("%u", cigar2->seq[c]);
    }
    printf("'\n");

    free(cigar_string);
    parasail_cigar_free(cigar);
    parasail_cigar_free(cigar2);

    parasail_result_free(result);

    /* let's see how the stats version of the same trace function
     * behaves */
    {
        char *dup = strdup_rpl(funcname);
        char *loc = strstr(dup, "trace");
        int score = 0;
        int end_query = 0;
        int end_ref = 0;
        int matches = 0;
        int similar = 0;
        int length = 0;
        if (NULL != loc) {
            /* this works because lengths are same */
            (void)strcpy(loc, "stats");
        }
        printf("stats function '%s'\n", dup);
        function = parasail_lookup_function(dup);
        if (NULL == function) {
            fprintf(stderr, "corresponding stats function not found.\n");
            exit(EXIT_FAILURE);
        }
        
        result = function(seqA, lena, seqB, lenb, open, extend, matrix);
        score = parasail_result_get_score(result);
        end_query = parasail_result_get_end_query(result);
        end_ref = parasail_result_get_end_ref(result);
        parasail_result_get_matches(result);
        parasail_result_get_similar(result);
        parasail_result_get_length(result);

        printf("Length:        %d\n", length);
        printf("Identity:   %d/%d\n", matches, length);
        printf("Similarity: %d/%d\n", similar, length);
        printf("Gaps:       %d/%d\n", -1, length);
        printf("Score:         %d\n", score);
        printf("end_query:     %d\n", end_query);
        printf("end_ref:       %d\n", end_ref);

        parasail_result_free(result);

        free(dup);
    }

    /* cleanup parsed sequences */
    parasail_sequences_free(sequences);

    return 0;
}

