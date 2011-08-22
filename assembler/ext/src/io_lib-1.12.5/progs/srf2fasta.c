/*
 * ======================================================================
 * This software has been created by Genome Research Limited (GRL).
 *
 * GRL hereby grants permission to use, copy, modify and distribute
 * this software and its documentation for non-commercial purposes
 * without fee at the user's own risk on the basis set out below.
 *
 * GRL neither undertakes nor accepts any duty whether contractual or
 * otherwise in connection with the software, its use or the use of
 * any derivative, and makes no representations or warranties, express
 * or implied, concerning the software, its suitability, fitness for
 * a particular purpose or non-infringement.
 *
 * In no event shall the authors of the software or GRL be responsible
 * or liable for any loss or damage whatsoever arising in any way
 * directly or indirectly out of the use of this software or its
 * derivatives, even if advised of the possibility of such damage.
 *
 * Our software can be freely distributed under the conditions set out
 * above, and must contain this copyright notice.
 * ======================================================================
 */

/*
 * This performs a linear (non-indexed) search for a trace in an SRF archive.
 *
 * It's not intended as a suitable production program or as a library of code
 * to use, but as a test and benchmark statistic.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>

#include <io_lib/Read.h>
#include <io_lib/misc.h>
#include <io_lib/ztr.h>
#include <io_lib/srf.h>

/* ------------------------------------------------------------------------ */

#define MAX_READ_LEN 10000
void ztr2fasta(ztr_t *z, char *name) {
    int i, nc;
    char buf[MAX_READ_LEN*2 + 512 + 6];
    char *seq = buf;
    ztr_chunk_t **chunks;

    /* Extract the sequence only */
    chunks = ztr_find_chunks(z, ZTR_TYPE_BASE, &nc);
    if (nc != 1) {
	fprintf(stderr, "Zero or greater than one BASE chunks found.\n");
	if (chunks)
	    free(chunks);
	return;
    }

    uncompress_chunk(z, chunks[0]);

    /* Construct fasta entry */
    *seq++ = '>';
    while (*name)
	*seq++ = *name++;
    *seq++ = '\n';

    for (i = 1; i < chunks[0]->dlength; i++) {
	char base = chunks[0]->data[i];
	if (base == '.')
	    *seq++ = 'N';
	else
	    *seq++ = base;
    }
    *seq++ = '\n';

    fwrite(buf, 1, seq - buf, stdout);
    free(chunks);

    return;
}

/* ------------------------------------------------------------------------ */
void usage(void) {
    fprintf(stderr, "Usage: srf2fasta [-C] archive_name\n");
    exit(1);
}

int main(int argc, char **argv) {
    char *ar_name;
    srf_t *srf;
    char name[512];
    ztr_t *ztr;
    int mask = 0, i;

    /* Parse args */
    for (i = 1; i < argc && argv[i][0] == '-'; i++) {
	if (!strcmp(argv[i], "-")) {
	    break;
	} else if (!strcmp(argv[i], "-C")) {
	    mask = SRF_READ_FLAG_BAD_MASK;
	} else {
	    usage();
	}
    }    

    if (i == argc) {
	usage();
    }
    ar_name = argv[i];

    if (NULL == (srf = srf_open(ar_name, "r"))) {
	perror(ar_name);
	return 4;
    }

    read_sections(READ_BASES);

#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif

    while (NULL != (ztr = srf_next_ztr(srf, name, mask))) {
	ztr2fasta(ztr, name);
	delete_ztr(ztr);
    }

    srf_destroy(srf, 1);

    return 0;
}
