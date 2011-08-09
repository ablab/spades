/*
 * Portions of this code have been derived from 454 Life Sciences Corporation's
 * getsff.c code (specifically the WriteSFFFile function.
 * It bears the following copyright notice:
 *
 * ------------------------------------------------------------
 * Copyright (c)[2001-2005] 454 Life Sciences Corporation. All Rights Reserved.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 *
 * IN NO EVENT SHALL LICENSOR BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided
 * that this copyright and notice appears in all copies.
 * ------------------------------------------------------------
 *
 * The remainder is Copyright Genome Research Limited (GRL) and is covered
 * by a BSD style license as described elsewhere in this source tree.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io_lib/Read.h"
#include "io_lib/xalloc.h"
#include "io_lib/sff.h"
#include "io_lib/misc.h"

/* -------------------------------------------------------------------------*/
/* General purpose decoding and mFILE reading functions */

/*
 * Unpacks the 31-byte fixed size part of the SFF common header.
 * It allocates memory for this and for the flow order and key, but does
 * not read the flow & key information (as this may not be in buf).
 * It also checks that the MAGIC and VERSION match as expected.
 *
 * Returns sff_common_header* on success
 *         NULL on failure
 */
sff_common_header *decode_sff_common_header(unsigned char *buf) {
    sff_common_header *h;

    if (NULL == (h = (sff_common_header *)xcalloc(1, sizeof(*h))))
	return NULL;

    h->magic           = be_int4(*(uint32_t *)(buf+0));
    memcpy(h->version, buf+4, 4);
    h->index_offset    = be_int8(*(uint64_t *)(buf+8));
    h->index_len       = be_int4(*(uint32_t *)(buf+16));
    h->nreads          = be_int4(*(uint32_t *)(buf+20));
    h->header_len      = be_int2(*(uint16_t *)(buf+24));
    h->key_len         = be_int2(*(uint16_t *)(buf+26));
    h->flow_len        = be_int2(*(uint16_t *)(buf+28));
    h->flowgram_format = be_int1(*(uint8_t  *)(buf+30));

    if (h->magic != SFF_MAGIC || memcmp(h->version, SFF_VERSION, 4)) {
	xfree(h);
	return NULL;
    }

    if (NULL == (h->flow = (char *)xmalloc(h->flow_len)))
	return free_sff_common_header(h), NULL;
    if (NULL == (h->key  = (char *)xmalloc(h->key_len)))
	return free_sff_common_header(h), NULL;

    return h;
}

/*
 * Encodes the data in 'h' to the file SFF representation. Buf should be
 * allocated to be 31 + h->flow_len + h->key_len + 8.
 *
 * Returns: the written length of buf
 */
int encode_sff_common_header(sff_common_header *h, unsigned char *buf) {
    int end;

    *(uint32_t *)(buf+0)  = be_int4(h->magic);
    memcpy(buf+4, h->version, 4);
    *(uint64_t *)(buf+8)  = be_int8(h->index_offset);
    *(uint32_t *)(buf+16) = be_int4(h->index_len);
    *(uint32_t *)(buf+20) = be_int4(h->nreads);
    *(uint16_t *)(buf+24) = be_int2(h->header_len);
    *(uint16_t *)(buf+26) = be_int2(h->key_len);
    *(uint16_t *)(buf+28) = be_int2(h->flow_len);
    *(uint8_t  *)(buf+30) = be_int1(h->flowgram_format);
    memcpy(buf+31, h->flow, h->flow_len);
    memcpy(buf+31+h->flow_len, h->key, h->key_len);
    end = 31+h->flow_len+h->key_len;
    memcpy(buf+end, "\0\0\0\0\0\0\0\0", ((end+7)&~7)-end);
    
    return (end+7)&~7;
}

/*
 * Reads a common header (including variable length components) from an mFILE.
 *
 * Returns the a pointer to the header on success
 *         NULL on failure
 */
sff_common_header *read_sff_common_header(mFILE *mf) {
    sff_common_header *h;
    unsigned char chdr[31];

    if (31 != mfread(chdr, 1, 31, mf))
	return NULL;
    h = decode_sff_common_header(chdr);

    if (h->flow_len != mfread(h->flow, 1, h->flow_len, mf))
	return free_sff_common_header(h), NULL;
    if (h->key_len != mfread(h->key , 1, h->key_len,  mf))
	return free_sff_common_header(h), NULL;

    /* Pad to 8 chars */
    mfseek(mf, (mftell(mf) + 7)& ~7, SEEK_SET);

    return h;
}

/*
 * Deallocates memory used by an SFF common header.
 */
void free_sff_common_header(sff_common_header *h) {
    if (!h)
	return;
    if (h->flow)
	xfree(h->flow);
    if (h->key)
	xfree(h->key);
    xfree(h);
}

/*
 * Unpacks the 16-byte fixed size part of the SFF read header.
 * It allocates memory for this and for the base calls, but does not
 * unpack these.
 *
 * Returns sff_read_header* on success
 *         NULL on failure
 */
sff_read_header *decode_sff_read_header(unsigned char *buf) {
    sff_read_header *h;

    if (NULL == (h = (sff_read_header *)xcalloc(1, sizeof(*h))))
	return NULL;

    h->header_len         = be_int2(*(uint16_t *)(buf+0));
    h->name_len           = be_int2(*(uint16_t *)(buf+2));
    h->nbases             = be_int4(*(uint32_t *)(buf+4));
    h->clip_qual_left     = be_int2(*(uint16_t *)(buf+8));
    h->clip_qual_right    = be_int2(*(uint16_t *)(buf+10));
    h->clip_adapter_left  = be_int2(*(uint16_t *)(buf+12));
    h->clip_adapter_right = be_int2(*(uint16_t *)(buf+14));

    if (NULL == (h->name  = (char *)xmalloc(h->name_len)))
	return free_sff_read_header(h), NULL;

    return h;
}

/*
 * Encodes the data in 'h' to the file SFF representation. Buf should be
 * allocated to be 16 + h->name_len + 8.
 *
 * Returns: the written length of buf
 */
int encode_sff_read_header(sff_read_header *h, unsigned char *buf) {
    int end;

    *(uint16_t *)(buf+0)  = be_int2(h->header_len);
    *(uint16_t *)(buf+2)  = be_int2(h->name_len);
    *(uint32_t *)(buf+4)  = be_int4(h->nbases);
    *(uint16_t *)(buf+8)  = be_int2(h->clip_qual_left);
    *(uint16_t *)(buf+10) = be_int2(h->clip_qual_right);
    *(uint16_t *)(buf+12) = be_int2(h->clip_adapter_left);
    *(uint16_t *)(buf+14) = be_int2(h->clip_adapter_right);
    memcpy(buf+16, h->name, h->name_len);
    end = 16+h->name_len;
    memcpy(buf+end, "\0\0\0\0\0\0\0\0", ((end+7)&~7)-end);
    
    return (end+7)&~7;
}

/*
 * Reads a read header (including variable length components) from an mFILE.
 *
 * Returns the a pointer to the header on success
 *         NULL on failure
 */
sff_read_header *read_sff_read_header(mFILE *mf) {
    sff_read_header *h;
    unsigned char rhdr[16];

    if (16 != mfread(rhdr, 1, 16, mf))
	return NULL;
    h = decode_sff_read_header(rhdr);

    if (h->name_len != mfread(h->name, 1, h->name_len, mf))
	return free_sff_read_header(h), NULL;
    
    /* Pad to 8 chars */
    mfseek(mf, (mftell(mf) + 7)& ~7, SEEK_SET);

    return h;
}

/*
 * Deallocates memory used by an SFF read header
 */
void free_sff_read_header(sff_read_header *h) {
    if (!h)
	return;
    if (h->name)
	xfree(h->name);
    free(h);
}

/*
 * Reads a read data block from an mFILE given a count of the number of
 * flows and basecalls (from the common header and read headers).
 *
 * Returns the a pointer to sff_read_data on success
 *         NULL on failure
 */
sff_read_data *read_sff_read_data(mFILE *mf, int nflows, int nbases) {
    sff_read_data *d;
    int i;

    if (NULL == (d = (sff_read_data *)xcalloc(1, sizeof(*d))))
	return NULL;

    if (NULL == (d->flowgram = (uint16_t *)xcalloc(nflows, 2)))
	return free_sff_read_data(d), NULL;
    if (nflows != mfread(d->flowgram, 2, nflows, mf))
	return free_sff_read_data(d), NULL;
    for (i = 0; i < nflows; i++)
	d->flowgram[i] = be_int2(d->flowgram[i]);

    if (NULL == (d->flow_index = (uint8_t *)xmalloc(nbases)))
	return free_sff_read_data(d), NULL;
    if (nbases != mfread(d->flow_index, 1, nbases, mf))
	return free_sff_read_data(d), NULL;

    if (NULL == (d->bases = (char *)xmalloc(nbases)))
	return free_sff_read_data(d), NULL;
    if (nbases != mfread(d->bases, 1, nbases, mf))
	return free_sff_read_data(d), NULL;

    if (NULL == (d->quality = (uint8_t *)xmalloc(nbases)))
	return free_sff_read_data(d), NULL;
    if (nbases != mfread(d->quality, 1, nbases, mf))
	return free_sff_read_data(d), NULL;

    /* Pad to 8 chars */
    mfseek(mf, (mftell(mf) + 7)& ~7, SEEK_SET);

    return d;
}

/*
 * Deallocates memory used by an SFF read data block
 */
void free_sff_read_data(sff_read_data *d) {
    if (!d)
	return;
    if (d->flowgram)
	xfree(d->flowgram);
    if (d->flow_index)
	xfree(d->flow_index);
    if (d->bases)
	xfree(d->bases);
    if (d->quality)
	xfree(d->quality);
    xfree(d);
}


/* -------------------------------------------------------------------------*/

/*
 * Reads an SFF file from an mFILE and decodes it to a Read struct.
 *
 * Returns Read* on success
 *         NULL on failure
 */
Read *mfread_sff(mFILE *mf) {
    int i, bpos;
    Read *r;
    sff_common_header *ch;
    sff_read_header *rh;
    sff_read_data *rd;

    /* Load the SFF contents */
    if (NULL == (ch = read_sff_common_header(mf)))
	return NULL;
    if (NULL == (rh = read_sff_read_header(mf))) {
	free_sff_common_header(ch);
	return NULL;
    }
    if (NULL == (rd = read_sff_read_data(mf, ch->flow_len, rh->nbases))) {
	free_sff_common_header(ch);
	free_sff_read_header(rh);
	return NULL;
    }

    /* Convert to Read struct */
    r = read_allocate(0,0);
    if (r->basePos) free(r->basePos);
    if (r->base)    free(r->base);
    if (r->prob_A)  free(r->prob_A);
    if (r->prob_C)  free(r->prob_C);
    if (r->prob_G)  free(r->prob_G);
    if (r->prob_T)  free(r->prob_T);

    r->nflows = ch->flow_len;
    r->flow_order = ch->flow; ch->flow = NULL;
    r->flow_raw = NULL;
    r->flow = (float *)malloc(r->nflows * sizeof(float));
    for (i = 0; i < r->nflows; i++) {
	r->flow[i] = rd->flowgram[i] / 100.0;
    }

    r->NBases = rh->nbases;
    r->basePos = (uint_2 *)calloc(r->NBases, 2);
    r->base    = rd->bases; rd->bases = NULL;
    r->prob_A  = (char *)calloc(r->NBases, 1);
    r->prob_C  = (char *)calloc(r->NBases, 1);
    r->prob_G  = (char *)calloc(r->NBases, 1);
    r->prob_T  = (char *)calloc(r->NBases, 1);

    bpos = 0;
    for (i=0; i < r->NBases; i++) {
	r->prob_A[i] = 0;
	r->prob_C[i] = 0;
	r->prob_G[i] = 0;
	r->prob_T[i] = 0;
	switch (r->base[i]) {
	case 'A':
	case 'a':
	    r->prob_A[i] = rd->quality[i];
	    break;
	case 'C':
	case 'c':
	    r->prob_C[i] = rd->quality[i];
	    break;
	case 'G':
	case 'g':
	    r->prob_G[i] = rd->quality[i];
	    break;
	case 'T':
	case 't':
	    r->prob_T[i] = rd->quality[i];
	    break;
	}

	bpos += rd->flow_index[i];
	r->basePos[i] = bpos;
    }

    r->leftCutoff = MAX(rh->clip_qual_left, rh->clip_adapter_left);
    r->rightCutoff = MIN(rh->clip_qual_right
			 ? rh->clip_qual_right
			 : r->NBases+1,
			 rh->clip_adapter_right
			 ? rh->clip_adapter_right
			 : r->NBases+1);

    free_sff_common_header(ch);
    free_sff_read_header(rh);
    free_sff_read_data(rd);

    return r;
}
