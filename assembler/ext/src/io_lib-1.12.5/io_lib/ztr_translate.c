#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <string.h>

#include "io_lib/ztr.h"
#include "io_lib/xalloc.h"
#include "io_lib/Read.h"

/* Return the A,C,G,T samples */
static char *ztr_encode_samples_4(ztr_t *z,
				  Read *r, int *nbytes, char **mdata,
				  int *mdbytes) {
    char *bytes;
    int i, j, t;

    if (!r->NPoints)
	return NULL;

    if ((z->header.version_major > 1 ||
	z->header.version_minor >= 2) && r->baseline) {
	/* 1.2 onwards */
	char buf[256];
	int blen;
	blen = sprintf(buf, "%d", r->baseline);
	*mdata = (char *)malloc(6+blen);
	*mdbytes = sprintf(*mdata, "OFFS%c%s", 0, buf) + 1;
    } else {
	*mdata = NULL;
	*mdbytes = 0;
    }

    bytes = (char *)xmalloc(r->NPoints * sizeof(TRACE)*4 + 2);
    for (i = 0, j = 2; i < r->NPoints; i++) {
	t = r->traceA[i];
	bytes[j++] = (t >> 8) & 0xff;
	bytes[j++] = (t >> 0) & 0xff;
    }
    for (i = 0; i < r->NPoints; i++) {
	t = r->traceC[i];
	bytes[j++] = (t >> 8) & 0xff;
	bytes[j++] = (t >> 0) & 0xff;
    }
    for (i = 0; i < r->NPoints; i++) {
	t = r->traceG[i];
	bytes[j++] = (t >> 8) & 0xff;
	bytes[j++] = (t >> 0) & 0xff;
    }
    for (i = 0; i < r->NPoints; i++) {
	t = r->traceT[i];
	bytes[j++] = (t >> 8) & 0xff;
	bytes[j++] = (t >> 0) & 0xff;
    }
    *nbytes = 4 * r->NPoints * sizeof(TRACE) + 2;

    bytes[0] = ZTR_FORM_RAW;
    bytes[1] = 0;
    return bytes;
}

static void ztr_decode_samples_4(ztr_t *z, ztr_chunk_t *chunk, Read *r) {
    int i, j;
    int maxTraceVal = 0;
    TRACE sample;
    unsigned char *bytes = (unsigned char *)chunk->data;
    int nbytes = chunk->dlength;

    bytes+=2;
    nbytes-=2;

    /* Store in the Read structure */
    r->NPoints = nbytes/8;
    if (r->traceA) xfree(r->traceA);
    if (r->traceC) xfree(r->traceC);
    if (r->traceG) xfree(r->traceG);
    if (r->traceT) xfree(r->traceT);
    r->traceA = (TRACE *)xmalloc(r->NPoints * sizeof(TRACE));
    r->traceC = (TRACE *)xmalloc(r->NPoints * sizeof(TRACE));
    r->traceG = (TRACE *)xmalloc(r->NPoints * sizeof(TRACE));
    r->traceT = (TRACE *)xmalloc(r->NPoints * sizeof(TRACE));

    for (i = j = 0; i < r->NPoints; i++, j+=2) {
	sample = (bytes[j] << 8) | bytes[j+1];
	r->traceA[i] = sample;
	if (maxTraceVal < sample)
	    maxTraceVal = sample;
    }
    for (i = 0; i < r->NPoints; i++, j+=2) {
	sample = (bytes[j] << 8) | bytes[j+1];
	r->traceC[i] = sample;
	if (maxTraceVal < sample)
	    maxTraceVal = sample;
    }
    for (i = 0; i < r->NPoints; i++, j+=2) {
	sample = (bytes[j] << 8) | bytes[j+1];
	r->traceG[i] = sample;
	if (maxTraceVal < sample)
	    maxTraceVal = sample;
    }
    for (i = 0; i < r->NPoints; i++, j+=2) {
	sample = (bytes[j] << 8) | bytes[j+1];
	r->traceT[i] = sample;
	if (maxTraceVal < sample)
	    maxTraceVal = sample;
    }

    r->maxTraceVal = maxTraceVal;
}

/* Return the [A,C,G,T] samples */
static char *ztr_encode_samples_common(ztr_t *z,
				       char ident[4], Read *r,
				       TRACE *data, int *nbytes,
				       char **mdata, int *mdbytes) {
    char *bytes;
    int i, j, t;

    if (!r->NPoints)
	return NULL;

    if (z->header.version_major > 1 ||
	z->header.version_minor >= 2) {
	/* 1.2 onwards */
	char buf[256];
	int blen;
	if (r->baseline) {
	    blen = sprintf(buf, "%d", r->baseline);
	    *mdata = (char *)malloc(16+blen);
	    *mdbytes = sprintf(*mdata, "TYPE%c%.*s%cOFFS%c%s",
			       0, 4, ident, 0,
			       0, buf) + 1;
	} else {
	    *mdata = (char *)malloc(10);
	    *mdbytes = sprintf(*mdata, "TYPE%c%.*s", 0, 4, ident) + 1;
	}
    } else {
	*mdata = (char *)malloc(4);
	*mdbytes = 4;
	(*mdata)[0] = ident[0];
	(*mdata)[1] = ident[1];
	(*mdata)[2] = ident[2];
	(*mdata)[3] = ident[3];
    }

    bytes = (char *)xmalloc(r->NPoints * sizeof(TRACE) + 2);
    for (i = 0, j = 2; i < r->NPoints; i++) {
	t = data[i];
	bytes[j++] = (t >> 8) & 0xff;
	bytes[j++] = (t >> 0) & 0xff;
    }
    *nbytes = r->NPoints * sizeof(TRACE) + 2;

    bytes[0] = ZTR_FORM_RAW;
    bytes[1] = 0;
    return bytes;
}


static char *ztr_encode_samples_A(ztr_t *z,
				  Read *r, int *nbytes, char **mdata,
				  int *mdbytes) {
    return ztr_encode_samples_common(z, "A\0\0", r, r->traceA,
				     nbytes, mdata, mdbytes);
}

static char *ztr_encode_samples_C(ztr_t *z,
				  Read *r, int *nbytes, char **mdata,
				  int *mdbytes) {
    return ztr_encode_samples_common(z, "C\0\0", r, r->traceC,
				     nbytes, mdata, mdbytes);
}

static char *ztr_encode_samples_G(ztr_t *z,
				  Read *r, int *nbytes, char **mdata,
				  int *mdbytes) {
    return ztr_encode_samples_common(z, "G\0\0", r, r->traceG,
				     nbytes, mdata, mdbytes);
}

static char *ztr_encode_samples_T(ztr_t *z,
				  Read *r, int *nbytes, char **mdata,
				  int *mdbytes) {
    return ztr_encode_samples_common(z, "T\0\0", r, r->traceT,
				     nbytes, mdata, mdbytes);
}

/* ARGSUSED */
static void ztr_decode_samples(ztr_t *z, ztr_chunk_t *chunk, Read *r) {
    int i, j;
    int maxTraceVal = 0;
    TRACE sample;
    unsigned char *bytes = (unsigned char *)chunk->data;
    int dlen = chunk->dlength;
    TRACE **lane, *lanex;
    char *type = ztr_lookup_mdata_value(z, chunk, "TYPE");
    
    if (!type)
	return;

    switch(type[0]) {
    case 'A':
	lane = &r->traceA;
	break;
    case 'C':
	lane = &r->traceC;
	break;
    case 'G':
	lane = &r->traceG;
	break;
    case 'T':
	lane = &r->traceT;
	break;
    default:
	return;
    }

    bytes+=2;
    dlen-=2;

    /* Store in the Read structure */
    r->NPoints = dlen/2;
    if (*lane)
	xfree(*lane);
    lanex = *lane = (TRACE *)xmalloc(r->NPoints * sizeof(TRACE));

    for (i = j = 0; i < r->NPoints; i++, j+=2) {
	sample = (bytes[j] << 8) | bytes[j+1];
	lanex[i] = sample;
	if (maxTraceVal < sample)
	    maxTraceVal = sample;
    }

    if (r->maxTraceVal < maxTraceVal)
	r->maxTraceVal = maxTraceVal;
}

/* Encode the the base calls */
static char *ztr_encode_bases(ztr_t *z,
			      Read *r, int *nbytes, char **mdata,
			      int *mdbytes) {
    char *bytes;

    if (!r->NBases)
	return NULL;
    
    *mdata = NULL;
    *mdbytes = 0;

    bytes = (char *)xmalloc(r->NBases + 1);
    memcpy(bytes+1, r->base, r->NBases);
    *nbytes = r->NBases+1;

    bytes[0] = ZTR_FORM_RAW;
    return bytes;
}

static void ztr_decode_bases(ztr_t *z, ztr_chunk_t *chunk, Read *r) {
    char *bytes = chunk->data;
    int nbytes = chunk->dlength;

    nbytes--;
    bytes++;

    r->NBases = nbytes;
    if (r->base) xfree(r->base);
    r->base = (char *)xmalloc(r->NBases+1);
    memcpy(r->base, bytes, r->NBases);
    r->base[r->NBases] = 0;

    /* Incase there isn't a clip chunk */
    r->leftCutoff = 0;
    r->rightCutoff = r->NBases+1;
}

/* Encode the base positions as 4 byte values */
static char *ztr_encode_positions(ztr_t *z,
				  Read *r, int *nbytes, char **mdata,
				  int *mdbytes) {
    char *bytes;
    int i, j;
    
    if ((!r->NPoints && !r->nflows) || !r->basePos || !r->NBases)
	return NULL;

    *mdata = NULL;
    *mdbytes = 0;

    bytes = (char *)xmalloc(r->NBases * 4 + 4);
    for (j = 4, i = 0; i < r->NBases; i++) {
	/*
	 * First 2 bytes are zero as currently r->basePos is 16-bit.
	 *
	 * bytes[j++] = (r->basePos[i] >> 24) & 0xff;
	 * bytes[j++] = (r->basePos[i] >> 16) & 0xff;
	 */
	bytes[j++] = 0;
	bytes[j++] = 0;
	bytes[j++] = (r->basePos[i] >>  8) & 0xff;
	bytes[j++] = (r->basePos[i] >>  0) & 0xff;
    }

    bytes[0] = ZTR_FORM_RAW;
    bytes[1] = 0; /* Dummy */
    bytes[2] = 0; /* Dummy */
    bytes[3] = 0; /* Dummy */

    *nbytes = j;
    return (char *)bytes;
}

static void ztr_decode_positions(ztr_t *z, ztr_chunk_t *chunk, Read *r) {
    int i, j;
    unsigned char *bytes = (unsigned char *)chunk->data;
    int nbytes = chunk->dlength;

    bytes+=4;
    nbytes-=4;

    r->NBases = nbytes/4;
    if (r->basePos) xfree(r->basePos);
    r->basePos = (uint_2 *)xmalloc(r->NBases * sizeof(*r->basePos));

    for (i = j = 0; j < nbytes; i++, j += 4) {
	r->basePos[i] =
	    (bytes[j+0] << 24) +
	    (bytes[j+1] << 16) +
	    (bytes[j+2] <<  8) +
	    (bytes[j+3] <<  0);
    }
}

/* Encode the main base confidence (called base) */
static char *ztr_encode_confidence_1(ztr_t *z,
				     Read *r, int *nbytes, char **mdata,
				     int *mdbytes)
{
    char *bytes;
    int i;

    /* Check that we have any confidence values first */
    if (!r->prob_A || !r->prob_C || !r->prob_G || !r->prob_T)
	return NULL;

    *mdata = NULL;
    *mdbytes = 0;

    /* Check that they're not all zero - will "normally" be quick */
    for (i = 0; i < r->NBases; i++) {
	if (r->prob_A[i]) break;
	if (r->prob_C[i]) break;
	if (r->prob_G[i]) break;
	if (r->prob_T[i]) break;
    }
    if (i == r->NBases)
	return NULL;

    /* Memory allocation */
    if (NULL == (bytes = xmalloc(r->NBases * sizeof(*bytes) + 1)))
	return NULL;

    /*
     * Encode probs for called bases.
     * Unknown base => average of prob_A, prob_C, prob_G and prob_T.
     */
    bytes++;
    for (i = 0; i < r->NBases; i++) {
	switch (r->base[i]) {
	case 'A':
	case 'a':
	    bytes[i] = r->prob_A[i];
	    break;
	case 'C':
	case 'c':
	    bytes[i] = r->prob_C[i];
	    break;
	case 'G':
	case 'g':
	    bytes[i] = r->prob_G[i];
	    break;
	case 'T':
	case 't':
	    bytes[i] = r->prob_T[i];
	    break;
	default:
	    bytes[i] = (r->prob_A[i] + r->prob_C[i] +
			r->prob_G[i] + r->prob_T[i]) / 4;
	    break;
	}
    }
    bytes--;

    *nbytes = r->NBases + 1;

    bytes[0] = ZTR_FORM_RAW;
    return bytes;
}

static int ztr_decode_confidence_1(ztr_t *z, ztr_chunk_t *chunk, Read *r) {
    char *bytes = chunk->data;
    int nbytes = chunk->dlength;
    int i;
    
    bytes++;
    nbytes--;

    /* Unpack confidence values; depends on base calls */
    if (!r->base)
	return -1;

    if (r->prob_A) xfree(r->prob_A);
    if (r->prob_C) xfree(r->prob_C);
    if (r->prob_G) xfree(r->prob_G);
    if (r->prob_T) xfree(r->prob_T);
    r->prob_A = (char *)xmalloc(r->NBases * sizeof(*r->prob_A));
    r->prob_C = (char *)xmalloc(r->NBases * sizeof(*r->prob_C));
    r->prob_G = (char *)xmalloc(r->NBases * sizeof(*r->prob_G));
    r->prob_T = (char *)xmalloc(r->NBases * sizeof(*r->prob_T));

    for (i = 0; i < r->NBases; i++) {
	switch (r->base[i]) {
	case 'A':
	case 'a':
	    r->prob_A[i] = bytes[i];
	    r->prob_C[i] = 0;
	    r->prob_G[i] = 0;
	    r->prob_T[i] = 0;
	    break;
	case 'C':
	case 'c':
	    r->prob_A[i] = 0;
	    r->prob_C[i] = bytes[i];
	    r->prob_G[i] = 0;
	    r->prob_T[i] = 0;
	    break;
	case 'G':
	case 'g':
	    r->prob_A[i] = 0;
	    r->prob_C[i] = 0;
	    r->prob_G[i] = bytes[i];
	    r->prob_T[i] = 0;
	    break;
	case 'T':
	case 't':
	    r->prob_A[i] = 0;
	    r->prob_C[i] = 0;
	    r->prob_G[i] = 0;
	    r->prob_T[i] = bytes[i];
	    break;
	default:
	    r->prob_A[i] = bytes[i];
	    r->prob_C[i] = bytes[i];
	    r->prob_G[i] = bytes[i];
	    r->prob_T[i] = bytes[i];
	}
    }

    return 0;
}

/* Encode the four main base confidences */
static char *ztr_encode_confidence_4(ztr_t *z,
				     Read *r, int *nbytes, char **mdata,
				     int *mdbytes)
{
    char *bytes;
    int i, j;

    /* Check that we have any confidence values first */
    if (!r->prob_A || !r->prob_C || !r->prob_G || !r->prob_T)
	return NULL;

    *mdata = NULL;
    *mdbytes = 0;

    /* Check that they're not all zero - will "normally" be quick */
    for (i = 0; i < r->NBases; i++) {
	if (r->prob_A[i]) break;
	if (r->prob_C[i]) break;
	if (r->prob_G[i]) break;
	if (r->prob_T[i]) break;
    }
    if (i == r->NBases)
	return NULL;

    /* Memory allocation */
    if (NULL == (bytes = xmalloc(4 * r->NBases * sizeof(*bytes) + 1)))
	return NULL;

    /*
     * Encode probs for called bases first
     * Unknown base = 'T'.
     */
    j = r->NBases;
    bytes++;
    for (i = 0; i < r->NBases; i++) {
	switch (r->base[i]) {
	case 'A':
	case 'a':
	    bytes[i  ] = r->prob_A[i];
	    bytes[j++] = r->prob_C[i];
	    bytes[j++] = r->prob_G[i];
	    bytes[j++] = r->prob_T[i];
	    break;
	case 'C':
	case 'c':
	    bytes[j++] = r->prob_A[i];
	    bytes[i  ] = r->prob_C[i];
	    bytes[j++] = r->prob_G[i];
	    bytes[j++] = r->prob_T[i];
	    break;
	case 'G':
	case 'g':
	    bytes[j++] = r->prob_A[i];
	    bytes[j++] = r->prob_C[i];
	    bytes[i  ] = r->prob_G[i];
	    bytes[j++] = r->prob_T[i];
	    break;
	default:
	    bytes[j++] = r->prob_A[i];
	    bytes[j++] = r->prob_C[i];
	    bytes[j++] = r->prob_G[i];
	    bytes[i  ] = r->prob_T[i];
	    break;
	}
    }
    bytes--;

    *nbytes = r->NBases * 4 + 1;

    bytes[0] = ZTR_FORM_RAW;
    return bytes;
}

static int ztr_decode_confidence_4(ztr_t *z, ztr_chunk_t *chunk, Read *r) {
    char *bytes = chunk->data;
    int nbytes = chunk->dlength;
    int i, j;
    
    bytes++;
    nbytes--;

    /* Unpack confidence values; depends on base calls */
    if (!r->base)
	return -1;

    if (r->prob_A) xfree(r->prob_A);
    if (r->prob_C) xfree(r->prob_C);
    if (r->prob_G) xfree(r->prob_G);
    if (r->prob_T) xfree(r->prob_T);
    r->prob_A = (char *)xmalloc(r->NBases * sizeof(*r->prob_A));
    r->prob_C = (char *)xmalloc(r->NBases * sizeof(*r->prob_C));
    r->prob_G = (char *)xmalloc(r->NBases * sizeof(*r->prob_G));
    r->prob_T = (char *)xmalloc(r->NBases * sizeof(*r->prob_T));

    j = r->NBases;
    for (i = 0; i < r->NBases; i++) {
	switch (r->base[i]) {
	case 'A':
	case 'a':
	    r->prob_A[i] = bytes[i];
	    r->prob_C[i] = bytes[j++];
	    r->prob_G[i] = bytes[j++];
	    r->prob_T[i] = bytes[j++];
	    break;
	case 'C':
	case 'c':
	    r->prob_A[i] = bytes[j++];
	    r->prob_C[i] = bytes[i];
	    r->prob_G[i] = bytes[j++];
	    r->prob_T[i] = bytes[j++];
	    break;
	case 'G':
	case 'g':
	    r->prob_A[i] = bytes[j++];
	    r->prob_C[i] = bytes[j++];
	    r->prob_G[i] = bytes[i];
	    r->prob_T[i] = bytes[j++];
	    break;
	default:
	    r->prob_A[i] = bytes[j++];
	    r->prob_C[i] = bytes[j++];
	    r->prob_G[i] = bytes[j++];
	    r->prob_T[i] = bytes[i];
	    break;
	}
    }

    return 0;
}

/* Encode the textual comments */
static char *ztr_encode_text(ztr_t *z,
			     Read *r, int *nbytes, char **mdata,
			     int *mdbytes) {
    char *bytes;
    int len, alen;
    int ident;
    int i, j;

    if (!r->info)
	return NULL;

    *mdata = NULL;
    *mdbytes = 0;

    /*
     * traditional Read comments are a single char * of ident=value lines.
     * The length of ident=valueXident=valueX (X = newline) if the same
     * as            ident0value0ident0value0 (0 = \0), although ztr has
     * a double \0 as terminator.
     */
    len = strlen(r->info);

    /* Allocate */
    alen = len + 3;
    bytes = xmalloc(alen);

    /* Copy */
    j = 0;
    bytes[j++] = ZTR_FORM_RAW;
    ident = 1;
    for (i = 0; i < len; i++) {
	switch (r->info[i]) {
	case '=':
	    if (ident) {
		ident = 0;
		bytes[j++] = 0;
	    } else {
		bytes[j++] = '=';
	    }
	    break;
	    
	case '\n':
	    if (ident) {
		/* Invalid Read info, but we'll carry on anyway. */
		if (j && bytes[j-1] != 0)
		    bytes[j++] = 0;
		else
		    break;
	    }
	    bytes[j++] = 0;
	    ident = 1;
	    break;

	default:
	    bytes[j++] = r->info[i];
	}

	if (j + 3 > alen) {
	    /* This can happen if we have Read idents without values */
	    alen += 100;
	    bytes = xrealloc(bytes, alen);
	}
    }
    if (j && bytes[j-1] != 0)
	bytes[j++] = 0; /* Must end in two nuls */
    bytes[j++] = 0;
    *nbytes = j;
    
    return bytes;
}

static void ztr_decode_text(Read *r, ztr_t *ztr) {
    int i;
    int nbytes = 0;
    char *iptr;

    /* Find length */
    for (i = 0; i < ztr->ntext_segments; i++) {
	if (ztr->text_segments[i].ident)
	    nbytes += strlen(ztr->text_segments[i].ident);
	if (ztr->text_segments[i].value)
	    nbytes += strlen(ztr->text_segments[i].value);
	nbytes += 2;
    }

    /* Allocate */
    if (r->info) xfree(r->info);
    r->info = (char *)xmalloc(nbytes+1);

    /* Convert */
    iptr = r->info;
    for (i = 0; i < ztr->ntext_segments; i++) {
	if (ztr->text_segments[i].ident && ztr->text_segments[i].value) {
	    int added = sprintf(iptr, "%s=%s\n", 
				ztr->text_segments[i].ident,
				ztr->text_segments[i].value);
	    iptr += added;
	}
    }
    *iptr = 0;
}

/* Encode the clip points */
static char *ztr_encode_clips(ztr_t *z,
			      Read *r, int *nbytes, char **mdata,
			      int *mdbytes) {
    char *bytes;

    if (!r->NBases)
	return NULL;

    if (r->leftCutoff == 0 && r->rightCutoff > r->NBases)
	return NULL;

    *mdata = NULL;
    *mdbytes = 0;

    /* Allocate */
    *nbytes = 9;
    bytes = xmalloc(9);

    /* Store */
    bytes[1] = (r->leftCutoff  >> 24) & 0xff;
    bytes[2] = (r->leftCutoff  >> 16) & 0xff;
    bytes[3] = (r->leftCutoff  >>  8) & 0xff;
    bytes[4] = (r->leftCutoff  >>  0) & 0xff;
    bytes[5] = (r->rightCutoff >> 24) & 0xff;
    bytes[6] = (r->rightCutoff >> 16) & 0xff;
    bytes[7] = (r->rightCutoff >>  8) & 0xff;
    bytes[8] = (r->rightCutoff >>  0) & 0xff;

    bytes[0] = ZTR_FORM_RAW;
    return bytes;
}

/* ARGSUSED */
static void ztr_decode_clips(ztr_t *z, ztr_chunk_t *chunk, Read *r) {
    char *bytes = chunk->data;

    r->leftCutoff =
	(((unsigned char)bytes[1]) << 24) +
	(((unsigned char)bytes[2]) << 16) +
	(((unsigned char)bytes[3]) <<  8) +
	(((unsigned char)bytes[4]) <<  0);

    r->rightCutoff =
	(((unsigned char)bytes[5]) << 24) +
	(((unsigned char)bytes[6]) << 16) +
	(((unsigned char)bytes[7]) <<  8) +
	(((unsigned char)bytes[8]) <<  0);
}

/* ARGSUSED */
static char *ztr_encode_flow_order(ztr_t *z,
				   Read *r, int *nbytes, char **mdata,
				   int *mdbytes) {
    char *bytes;

    if (!r->flow_order || !r->nflows)
	return NULL;

    bytes = (char *)xmalloc(r->nflows+1);
    *nbytes = r->nflows+1;
    bytes[0] = ZTR_FORM_RAW;

    memcpy(bytes+1, r->flow_order, r->nflows);

    return bytes;
}

static void ztr_decode_flow_order(ztr_t *z, ztr_chunk_t *chunk, Read *r) {
    char *bytes = chunk->data;
    int nbytes = chunk->dlength;

    nbytes--;
    bytes++;

    r->nflows = nbytes;
    r->flow_order = (char *)xmalloc(r->nflows+1);
    memcpy(r->flow_order, bytes, r->nflows);
    r->flow_order[r->nflows] = 0;
}

static char *ztr_encode_flow_proc(ztr_t *z,
				  Read *r, int *nbytes, char **mdata,
				  int *mdbytes) {
    char *bytes;
    int i, j;
    float *data;

    if (!r->flow_order || !r->nflows)
	return NULL;

    data = r->flow;

    /* Meta-data */
    if (z->header.version_major > 1 ||
	z->header.version_minor >= 2) {
	/* 1.2 onwards */
	*mdata = (char *)malloc(10);
	*mdbytes = 10;
	sprintf(*mdata, "TYPE%cPYNO", 0);
    } else {
	*mdata = (char *)malloc(4);
	*mdbytes = 4;
	(*mdata)[0] = 'P';
	(*mdata)[1] = 'Y';
	(*mdata)[2] = 'N';
	(*mdata)[3] = 'O';
    }

    /* floats themselves, scaled */
    bytes = (char *)xmalloc(r->nflows*2+2);
    *nbytes = r->nflows*2+2;
    bytes[0] = ZTR_FORM_RAW;
    bytes[1] = 0;

    for (i = 0, j = 2; i < r->nflows; i++) {
	signed int t = data[i] * 100 + 0.49999;
	bytes[j++] = (t >> 8) & 0xff;
	bytes[j++] = (t >> 0) & 0xff;
    }

    return bytes;
}

/* ARGSUSED */
static void ztr_decode_flow_proc(ztr_t *z, ztr_chunk_t *chunk, Read *r) {
    int i, j;
    unsigned char *bytes = (unsigned char *)chunk->data;
    int dlen = chunk->dlength;

    bytes+=2;
    dlen-=2;

    /* Store in the Read structure */
    r->nflows = dlen/2;
    r->flow = (float *)xcalloc(r->nflows, sizeof(float));
    for (i = j = 0; i < r->nflows; i++, j+=2) {
	float sample = ((bytes[j] << 8) | bytes[j+1]) / 100.0;
	r->flow[i] = sample;
    }
}

static char *ztr_encode_flow_raw(ztr_t *z,
				 Read *r, int *nbytes, char **mdata,
				 int *mdbytes) {
    char *bytes;
    int i, j;
    unsigned int *data;

    if (!r->flow_raw || !r->nflows)
	return NULL;

    data = r->flow_raw;

    /* Meta-data */
    if (z->header.version_major > 1 ||
	z->header.version_minor >= 2) {
	/* 1.2 onwards */
	*mdata = (char *)malloc(10);
	*mdbytes = 10;
	sprintf(*mdata, "TYPE%cPYRW", 0);
    } else {
	*mdata = (char *)malloc(4);
	*mdbytes = 4;
	(*mdata)[0] = 'P';
	(*mdata)[1] = 'Y';
	(*mdata)[2] = 'R';
	(*mdata)[3] = 'W';
    }

    /* floats themselves, scaled */
    bytes = (char *)xmalloc(r->nflows*2+2);
    *nbytes = r->nflows*2+2;
    bytes[0] = ZTR_FORM_RAW;
    bytes[1] = 0;

    for (i = 0, j = 2; i < r->nflows; i++) {
	int t = data[i];
	bytes[j++] = (t >> 8) & 0xff;
	bytes[j++] = (t >> 0) & 0xff;
    }

    return bytes;
}

/* ARGSUSED */
static void ztr_decode_flow_raw(ztr_t *z, ztr_chunk_t *chunk, Read *r) {
    int i, j;
    unsigned char *bytes = (unsigned char *)chunk->data;
    int dlen = chunk->dlength;

    bytes+=2;
    dlen-=2;

    /* Store in the Read structure */
    r->nflows = dlen/2;
    r->flow_raw = (unsigned int *)xcalloc(r->nflows, sizeof(*r->flow_raw));
    for (i = j = 0; i < r->nflows; i++, j+=2) {
	unsigned int sample = (bytes[j] << 8) | bytes[j+1];
	r->flow_raw[i] = sample;
    }
}

/*
 * read2ztr
 *
 * Converts an io_lib "Read" structure to a ztr_t structure.
 *
 * Arguments:
 *	r		A pointer to the "Read" structure to convert from
 *
 * Returns:
 *	Success:  A pointer to the ztr_t struct.
 *	Failure:  NULL
 */
ztr_t *read2ztr(Read *r) {
    ztr_t *ztr;
    int i, j, nbytes, mdbytes;
    char *bytes;
    char *mdata;

#define DO_SMP4
    int chunk_type[] = {
#ifdef DO_SMP4
	ZTR_TYPE_SMP4,
#else
	ZTR_TYPE_SAMP,
	ZTR_TYPE_SAMP,
	ZTR_TYPE_SAMP,
	ZTR_TYPE_SAMP,
#endif
	ZTR_TYPE_BASE,
	ZTR_TYPE_BPOS,
	ZTR_TYPE_CNF4,
	ZTR_TYPE_TEXT,
	ZTR_TYPE_CLIP,
	
	ZTR_TYPE_FLWO,
	ZTR_TYPE_SAMP,
	ZTR_TYPE_SAMP,
    };

    char *(*chunk_func[])(ztr_t *z, Read *r, int *nbytes, char **mdata, int *mdbytes) = {
#ifdef DO_SMP4
	ztr_encode_samples_4,
#else
	ztr_encode_samples_A,
	ztr_encode_samples_C,
	ztr_encode_samples_G,
	ztr_encode_samples_T,
#endif
	ztr_encode_bases,
	ztr_encode_positions,
	ztr_encode_confidence_4,
	ztr_encode_text,
	ztr_encode_clips,

	ztr_encode_flow_order,
	ztr_encode_flow_proc,
	ztr_encode_flow_raw,
    };

    if (NULL == (ztr = new_ztr()))
	return NULL;

    /* Create a header record */
    memcpy(ztr->header.magic, ZTR_MAGIC, 8);
    ztr->header.version_major = ZTR_VERSION_MAJOR;
    ztr->header.version_minor = ZTR_VERSION_MINOR;

    /* Alloc chunks (max number) */
    ztr->nchunks = sizeof(chunk_type)/sizeof(*chunk_type);
    ztr->chunk = (ztr_chunk_t *)xmalloc(ztr->nchunks *
					sizeof(ztr_chunk_t));
    if (NULL == ztr->chunk)
	return NULL;

    /* Create the chunks */
    for (j = i = 0; i < ztr->nchunks; i++) {
	/* char str[5]; */

	bytes = chunk_func[i](ztr, r, &nbytes, &mdata, &mdbytes);
	if (!bytes)
	    continue;

	/*
	fprintf(stderr, "block %.4s length %d\n",
		ZTR_BE2STR(chunk_type[i], str), nbytes);
	*/
	ztr->chunk[j].type     = chunk_type[i];
	ztr->chunk[j].mdlength = mdbytes;
	ztr->chunk[j].mdata    = mdata;
	ztr->chunk[j].dlength  = nbytes;
	ztr->chunk[j].data     = bytes;
	ztr->chunk[j].ztr_owns = 1;

	j++;
    }
    ztr->nchunks = j;

    /*
     * Experiments show that typically a double delta does
     * better than a single delta for 8-bit data, and the other
     * way around for 16-bit data
     */
    ztr->delta_level = r->maxTraceVal < 256 ? 2 : 3;
    
    return ztr;
}

/*
 * ztr2read
 *
 * Converts an ztr_t structure to an io_lib "Read" structure.
 *
 * Arguments:
 *	ztr		A pointer to the ztr structure to convert from
 *
 * Returns:
 *	Success:  A pointer to the Read struct.
 *	Failure:  NULL
 */
Read *ztr2read(ztr_t *ztr) {
    Read *r;
    int i;
    int done_conf = 0, done_pos = 0;
    int sections = read_sections(0);

    /* Allocate */
    r = read_allocate(0, 0);

    if (NULLRead == r)
	return NULLRead;

    /* Proces text chunks - makes conversion easier */
    if (sections & READ_COMMENTS) {
	ztr_process_text(ztr);
	ztr_decode_text(r, ztr);
    }

    /* Iterate around each known chunk type turning into the Read elements */
    for (i = 0; i < ztr->nchunks; i++) {
	switch (ztr->chunk[i].type) {
	case ZTR_TYPE_SMP4:
	    if (sections & READ_SAMPLES) {
		char *offs = ztr_lookup_mdata_value(ztr, &ztr->chunk[i], "OFFS");
		char *type = ztr_lookup_mdata_value(ztr, &ztr->chunk[i], "TYPE");
		if (!type || 0 == strcmp(type, "PROC")) {
		//if (type && 0 == strcmp(type, "SLXI")) {
		    uncompress_chunk(ztr, &ztr->chunk[i]);
		    ztr_decode_samples_4(ztr, &ztr->chunk[i], r);

		    if (offs)
			r->baseline = atoi(offs);
		}
	    }
	    break;

	case ZTR_TYPE_SAMP: 
	    if (sections & READ_SAMPLES) {
		char *type = ztr_lookup_mdata_value(ztr, &ztr->chunk[i], "TYPE");
		char *offs = ztr_lookup_mdata_value(ztr, &ztr->chunk[i], "OFFS");
		uncompress_chunk(ztr, &ztr->chunk[i]);
		if (type && 0 == strcmp(type, "PYRW"))
		    ztr_decode_flow_raw(ztr, &ztr->chunk[i], r);
		else if (type && 0 == strcmp(type, "PYNO"))
		    ztr_decode_flow_proc(ztr, &ztr->chunk[i], r);
		else if (type &&
			 (0 == strcmp(type, "A") ||
			  0 == strcmp(type, "C") ||
			  0 == strcmp(type, "G") ||
			  0 == strcmp(type, "T"))) {
		    ztr_decode_samples(ztr, &ztr->chunk[i], r);

		    if (offs)
			r->baseline = atoi(offs);
		}
	    }
	    break;

	case ZTR_TYPE_BASE:
	    if (sections & READ_BASES) {
		uncompress_chunk(ztr, &ztr->chunk[i]);
		ztr_decode_bases(ztr, &ztr->chunk[i], r);
	    }
	    break;

	case ZTR_TYPE_BPOS:
	    if (sections & READ_BASES) {
		uncompress_chunk(ztr, &ztr->chunk[i]);
		ztr_decode_positions(ztr, &ztr->chunk[i], r);
		done_pos++;
	    }
	    break;

	case ZTR_TYPE_CNF4:
	    if (sections & READ_BASES) {
		uncompress_chunk(ztr, &ztr->chunk[i]);
		ztr_decode_confidence_4(ztr, &ztr->chunk[i], r);
		done_conf++;
	    }
	    break;

	case ZTR_TYPE_CNF1:
	    if (sections & READ_BASES) {
		uncompress_chunk(ztr, &ztr->chunk[i]);
		ztr_decode_confidence_1(ztr, &ztr->chunk[i], r);
		done_conf++;
	    }
	    break;

	case ZTR_TYPE_TEXT:
	    /* Skip - already did this; see ztr_process_text */
	    break;

	case ZTR_TYPE_CLIP:
	    if (sections & READ_BASES) {
		uncompress_chunk(ztr, &ztr->chunk[i]);
		ztr_decode_clips(ztr, &ztr->chunk[i], r);
	    }
	    break;

	case ZTR_TYPE_FLWO:
	    if (sections & READ_SAMPLES) {
		uncompress_chunk(ztr, &ztr->chunk[i]);
		ztr_decode_flow_order(ztr, &ztr->chunk[i], r);
	    }
	    break;
	}
    }

    /* Handle the case when we have no confidence values */
    if (!done_conf && r->NBases > 0) {
	r->prob_A = (char *)xrealloc(r->prob_A, r->NBases);
	r->prob_C = (char *)xrealloc(r->prob_C, r->NBases);
	r->prob_G = (char *)xrealloc(r->prob_G, r->NBases);
	r->prob_T = (char *)xrealloc(r->prob_T, r->NBases);
	memset(r->prob_A, 0, r->NBases);
	memset(r->prob_C, 0, r->NBases);
	memset(r->prob_G, 0, r->NBases);
	memset(r->prob_T, 0, r->NBases);
    }

    /* Handle the case when we have no BPOS chunk */
    if (!done_pos && r->NBases > 0) {
	r->basePos = (uint_2 *)xrealloc(r->basePos, r->NBases * 2);
	for (i = 0; i < r->NBases; i++)
	    r->basePos[i] = i;
    }

    r->format = TT_ZTR;

    return r;
}
