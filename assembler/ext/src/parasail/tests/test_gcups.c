#include "config.h"

/* getopt needs _POSIX_C_SOURCE 2 */
#define _POSIX_C_SOURCE 2

#include <ctype.h>
#include <limits.h>
#include <math.h>
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

#include "parasail/parasail/io.h"


static inline unsigned long binomial_coefficient(unsigned long n, unsigned long k)
{
    /* from http://blog.plover.com/math/choose.html */
    unsigned long r = 1;
    unsigned long d;
    if (k > n) {
        return 0;
    }
    for (d = 1; d <= k; d++) {
        r *= n--;
        r /= d;
    }
    return r;
}

static inline void k_combination2(unsigned long pos, unsigned long *a, unsigned long *b)
{
    double s;
    double i = floor(sqrt(2.0 * pos)) - 1.0;
    if (i <= 1.0) {
        i = 1.0;
    }
    s = i * (i - 1.0) / 2.0;
    while (pos - s >= i) {
        s += i;
        i += 1;
    }
    *a = (unsigned long)(pos - s);
    *b = (unsigned long)(i);
}

int main(int argc, char **argv)
{
    unsigned long i = 0;
    size_t seq_count = 10;
    size_t limit = 0;
    parasail_sequences_t *sequences = NULL;
    char *filename = NULL;
    int c = 0;
    int distribution = 0;

    while ((c = getopt(argc, argv, "f:d")) != -1) {
        switch (c) {
            case 'd':
                distribution = 1;
                break;
            case 'f':
                filename = optarg;
                break;
            case '?':
                if (optopt == 'f' || optopt == 'n') {
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
        fprintf(stderr, "missing filename\n");
        exit(1);
    }

    if (distribution) {
        for (i=0; i<seq_count; ++i) {
            printf("%lu\n", (unsigned long)sequences->seqs[i].seq.l);
        }
    }
    else {
        limit = binomial_coefficient(seq_count, 2);
        printf("%lu choose 2 is %lu\n",
                (unsigned long)seq_count, (unsigned long)limit);

        {
            unsigned long a=0;
            unsigned long b=1;
            unsigned long work=0;
            unsigned long columns=0;
            for (i=0; i<limit; ++i) {
                k_combination2(i, &a, &b);
                work += sequences->seqs[a].seq.l*sequences->seqs[b].seq.l;
                columns += sequences->seqs[b].seq.l;
            }
            printf("work=%lu\n", work);
            printf("columns=%lu\n", columns);
        }
    }

    return 0;
}

