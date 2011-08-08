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
#include <string.h>
#include <fcntl.h>
#include <io_lib/Read.h>
#include <io_lib/misc.h>
#include <io_lib/ztr.h>
#include <io_lib/srf.h>

/* ------------------------------------------------------------------------ */
/*
 * Looks for a trace name in an SRF archive and returns the binary contents
 * if found, or NULL if not.
 */
mFILE *find_reading(srf_t *srf, char *tr_name) {
    do {
	int type;

	switch(type = srf_next_block_type(srf)) {
	case -1:
	    /* EOF */
	    return NULL;

	case SRFB_CONTAINER:
	    if (0 != srf_read_cont_hdr(srf, &srf->ch))
		return NULL;
	    break;

	case SRFB_XML:
	    if (0 != srf_read_xml(srf, &srf->xml))
		return NULL;
	    break;
	    

	case SRFB_TRACE_HEADER: {
	    /* off_t pos = ftell(srf->fp); */

	    if (0 != srf_read_trace_hdr(srf, &srf->th))
		return NULL;

#if 0
	    /*
	     * If the name prefix doesn't match tr_name then skip this entire
	     * block.
	     */
	    if (0 != strncmp(tr_name, srf->th.id_prefix,
			     strlen(srf->th.id_prefix)) &&
		0 != srf->th.next_block_offset) {
		fseek(srf->fp, pos + srf->th.next_block_offset, SEEK_SET);
	    }
#endif
	    break;
	}

	case SRFB_TRACE_BODY: {
	    mFILE *mf = mfcreate(NULL, 0);
	    srf_trace_body_t tb;
	    char name[512];

	    if (!mf || 0 != srf_read_trace_body(srf, &tb, 0))
		return NULL;

	    sprintf(name, "%s%s", srf->th.id_prefix, tb.read_id);
	    if (strcmp(name, tr_name)) {
		mfdestroy(mf);
		if (tb.trace)
		    free(tb.trace);
		continue;
	    }

	    if (srf->th.trace_hdr_size)
		mfwrite(srf->th.trace_hdr, 1,
			srf->th.trace_hdr_size, mf);
	    if (tb.trace_size)
		mfwrite(tb.trace, 1, tb.trace_size, mf);
	    if (tb.trace)
		free(tb.trace);
	    mrewind(mf);
	    return mf;
	}

	case SRFB_INDEX: {
	    off_t pos = ftello(srf->fp);
	    srf_read_index_hdr(srf, &srf->hdr, 1);

	    /* Skip the index body */
	    fseeko(srf->fp, pos + srf->hdr.size, SEEK_SET);
	    break;
	}

 	case SRFB_NULL_INDEX:
 	    break;

	default:
	    fprintf(stderr, "Block of unknown type '%c'. Aborting\n", type);
	    return NULL;
	}
    } while (1);

    return NULL;
}

/* ------------------------------------------------------------------------ */
int main(int argc, char **argv) {
    char *ar_name, *tr_name;
    mFILE *mf;
    srf_t *srf;

    if (argc != 3) {
	fprintf(stderr, "Usage: extract_linear_srf archive_name trace_name\n");
	return 1;
    }
    ar_name = argv[1];
    tr_name = argv[2];

    if (NULL == (srf = srf_open(ar_name, "r"))) {
	perror(ar_name);
	return 4;
    }

    if (NULL == (mf = find_reading(srf, tr_name))) {
	fprintf(stderr, "%s not found in archive\n", tr_name);
	return 3;
    }

#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif

    fwrite(mf->data, 1, mf->size, stdout);
    return 0;
}
