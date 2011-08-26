/*
 * Copyright (c) Medical Research Council 1994. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written by James Bonfield, Simon Dear, Rodger Staden,
 * as part of the Staden Package at the MRC Laboratory of Molecular
 * Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 */

/* 
    Title:       write_scf.c

    Purpose:	 Output of Standard Chromatogram Format sequences
    Last update: August 18 1994

    Change log:
    4 Feb 1992, Now draft proposal version 2
    23 Nov 92,  SCF 2.0 + LaDeana's changes
    11 Aug 93, Version 2.01 containing confidence values
    18 Aug 1994  Renamed from  writeSCF.c; now purely SCF IO (no Seq structs)

    Oct 95 major rewrite to make files more easily compressed.
    gzip now gets files to around 40% of original
    Version raised to 3.00
     * We store in order:
     *     Header
     *     Samples
     *     Bases
     *     Comments
     *     Private

     Two main types of change: 
     1: write data in lane order instead of all lanes together
     eg write Sample values for A, then Sample values for C, etc. 

     2: where appropriate write delta delta values instead of complete ones.
     ie write the differences in the differences between successive values

*/


static int scf_version = 3;

/* ---- Imports ---- */


#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <ctype.h>
#include <stdio.h>    /* IMPORT: fopen, fclose, fseek, ftell, fgetc,
		                 EOF */
#include <string.h>
#include "io_lib/scf.h"      /* IMPORT: scf structures */
#include "io_lib/mach-io.h"  /* IMPORT: be_write_int_1, be_write_int_2, be_write_int_4 */
#include "io_lib/xalloc.h"

#include "io_lib/stdio_hack.h"

/* ---- Exports ---- */


int write_scf_header(FILE *fp, Header *h)
{
    int i;

    if (be_write_int_4(fp,&h->magic_number)==0)         return -1;
    if (be_write_int_4(fp,&h->samples)==0)              return -1;
    if (be_write_int_4(fp,&h->samples_offset)==0)       return -1;
    if (be_write_int_4(fp,&h->bases)==0)                return -1;
    if (be_write_int_4(fp,&h->bases_left_clip)==0)      return -1;
    if (be_write_int_4(fp,&h->bases_right_clip)==0)     return -1;
    if (be_write_int_4(fp,&h->bases_offset)==0)         return -1;
    if (be_write_int_4(fp,&h->comments_size)==0)        return -1;
    if (be_write_int_4(fp,&h->comments_offset)==0)      return -1;
    if (fwrite(h->version,sizeof(h->version),1,fp)!=1)  return -1;
    if (be_write_int_4(fp,&h->sample_size)==0)          return -1;
    if (be_write_int_4(fp,&h->code_set)==0)             return -1;
    if (be_write_int_4(fp,&h->private_size)==0)         return -1;
    if (be_write_int_4(fp,&h->private_offset)==0)       return -1;
    for (i=0;i<18;i++)
	if (be_write_int_4(fp,&h->spare[i])==0)         return -1;

    return 0;
}


int write_scf_sample1(FILE *fp, Samples1 *s)
{
    uint_1 buf[4];

    buf[0] = s->sample_A;
    buf[1] = s->sample_C;
    buf[2] = s->sample_G;
    buf[3] = s->sample_T;
    if (4 != fwrite(buf, 1, 4, fp)) return -1;

    return 0;
}


int write_scf_sample2(FILE *fp, Samples2 *s)
{
    uint_2 buf[4];

    buf[0] = be_int2(s->sample_A);
    buf[1] = be_int2(s->sample_C);
    buf[2] = be_int2(s->sample_G);
    buf[3] = be_int2(s->sample_T);
    if (4 != fwrite(buf, 2, 4, fp)) return -1;

    return 0;
}


int write_scf_samples1(FILE *fp, Samples1 *s, size_t num_samples) {
    size_t i;

    for (i = 0; i < num_samples; i++) {
	if (-1 == write_scf_sample1(fp, &(s[i])))
	    return -1;
    }

    return 0;
}


int write_scf_samples2(FILE *fp, Samples2 *s, size_t num_samples) {
    size_t i;

    for (i = 0; i < num_samples; i++) {
	if (-1 == write_scf_sample2(fp, &(s[i])))
	    return -1;
    }

    return 0;
}


int write_scf_samples31(FILE *fp, Samples1 *s, size_t num_samples) {
    size_t i;
    int1 *samples_out;

    if (!num_samples)
	return 0;

    if ( ! (samples_out = (int1 *)xmalloc(num_samples * 
					    sizeof(int1)))) {
	return -1;
    }

    for (i = 0; i < num_samples; i++) {
	samples_out[i] = (&s[i])->sample_A;
    }
    scf_delta_samples1 ( samples_out, num_samples, 1);
    if (num_samples != fwrite(samples_out, 1, num_samples, fp)) {
	xfree(samples_out);
	return -1;
    }

    for (i = 0; i < num_samples; i++) {
	samples_out[i] = (&s[i])->sample_C;
    }
    scf_delta_samples1 ( samples_out, num_samples, 1);
    if (num_samples != fwrite(samples_out, 1, num_samples, fp)) {
	xfree(samples_out);
	return -1;
    }

    for (i = 0; i < num_samples; i++) {
	samples_out[i] = (&s[i])->sample_G;
    }
    scf_delta_samples1 ( samples_out, num_samples, 1);
    if (num_samples != fwrite(samples_out, 1, num_samples, fp)) {
	xfree(samples_out);
	return -1;
    }

    for (i = 0; i < num_samples; i++) {
	samples_out[i] = (&s[i])->sample_T;
    }
    scf_delta_samples1 ( samples_out, num_samples, 1);
    if (num_samples != fwrite(samples_out, 1, num_samples, fp)) {
	xfree(samples_out);
	return -1;
    }

    xfree(samples_out);
    return 0;
}

int write_scf_samples32(FILE *fp, Samples2 *s, size_t num_samples) {
    size_t i;
    uint2 *samples_out;

    if (!num_samples)
	return 0;

    if ( ! (samples_out = (uint2 *)xmalloc(num_samples * sizeof(uint2)))) {
	return -1;
    }


    for (i = 0; i < num_samples; i++) {
	samples_out[i] = (&s[i])->sample_A;
    }
    scf_delta_samples2 ( samples_out, num_samples, 1);
#ifdef SP_LITTLE_ENDIAN
    for (i = 0; i < num_samples; i++) {
	samples_out[i] = be_int2(samples_out[i]);
    }
#endif
    if (num_samples != fwrite(samples_out, 2, num_samples, fp)) return -1;


    for (i = 0; i < num_samples; i++) {
	samples_out[i] = (&s[i])->sample_C;
    }
    scf_delta_samples2 ( samples_out, num_samples, 1);
#ifdef SP_LITTLE_ENDIAN
    for (i = 0; i < num_samples; i++) {
	samples_out[i] = be_int2(samples_out[i]);
    }
#endif
    if (num_samples != fwrite(samples_out, 2, num_samples, fp)) return -1;


    for (i = 0; i < num_samples; i++) {
	samples_out[i] = (&s[i])->sample_G;
    }
    scf_delta_samples2 ( samples_out, num_samples, 1);
#ifdef SP_LITTLE_ENDIAN
    for (i = 0; i < num_samples; i++) {
	samples_out[i] = be_int2(samples_out[i]);
    }
#endif
    if (num_samples != fwrite(samples_out, 2, num_samples, fp)) return -1;


    for (i = 0; i < num_samples; i++) {
	samples_out[i] = (&s[i])->sample_T;
    }
    scf_delta_samples2 ( samples_out, num_samples, 1);
#ifdef SP_LITTLE_ENDIAN
    for (i = 0; i < num_samples; i++) {
	samples_out[i] = be_int2(samples_out[i]);
    }
#endif
    if (num_samples != fwrite(samples_out, 2, num_samples, fp)) return -1;


    xfree(samples_out);
    return 0;
}


int write_scf_base(FILE *fp, Bases *b)
{
    uint_1 buf[12];

    ((uint_4 *)buf)[0] = be_int4(b->peak_index);
    buf[4] = b->prob_A;
    buf[5] = b->prob_C;
    buf[6] = b->prob_G;
    buf[7] = b->prob_T;
    buf[8] = b->base;
    buf[9] = b->spare[0];
    buf[10] = b->spare[1];
    buf[11] = b->spare[2];

    if (12 != fwrite(buf, 1, 12, fp)) return -1;

    return 0;
}


int write_scf_bases(FILE *fp, Bases *b, size_t num_bases)
{
    size_t i;

    for (i = 0; i < num_bases; i++) {
	if (-1 == write_scf_base(fp, &(b[i])))
	    return -1;
    }

    return 0;
}

int write_scf_bases3(FILE *fp, Bases *b, size_t num_bases)
{
    size_t i;
    uint_4 *buf4;
    uint_1 *buf1;

    if (NULL == (buf4 = (uint_4 *)xmalloc(1 + 4 * num_bases)))
	return -1;

    if (NULL == (buf1 = (uint_1 *)xmalloc(1 + 8 * num_bases))) {
	xfree(buf4);
	return -1;
    }

    for (i = 0; i < num_bases; i++) {
	buf4[i] = be_int4((&b[i])->peak_index);
    }
    fwrite(buf4, 4, num_bases, fp);
    
    for (i=0; i < num_bases; i++) {
	buf1[i            ] = (&b[i])->prob_A;
	buf1[i+  num_bases] = (&b[i])->prob_C;
	buf1[i+2*num_bases] = (&b[i])->prob_G;
	buf1[i+3*num_bases] = (&b[i])->prob_T;
	buf1[i+4*num_bases] = (&b[i])->base;
	buf1[i+5*num_bases] = (&b[i])->spare[0];
	buf1[i+6*num_bases] = (&b[i])->spare[1];
	buf1[i+7*num_bases] = (&b[i])->spare[2];
    }
    if (8 * num_bases != (fwrite(buf1, 1, 8 * num_bases, fp))) {
	xfree(buf1);
	xfree(buf4);
	return -1;
    }

    xfree(buf1);
    xfree(buf4);
    return 0;
}


int write_scf_comment(FILE *fp, Comments *c, size_t s)
{
    if (fwrite(c, 1, s, fp) != s) return -1;

    return 0;
}



/*
 * Request which (major) version of scf to use when writing.
 * Defaults to the latest. Currently suitable fields are
 * 2 and 3.
 *
 * Returns 0 for success, -1 for failure.
 */
int set_scf_version(int version) {
    if (version != 2 && version != 3)
	return -1;

    scf_version = version;
    return 0;
}

/*
 * Write Seq out as a .scf file to the 'fp' FILE *
 */
int fwrite_scf(Scf *scf, FILE *fp) {
    uint_4 size;
    int err;

    /*
     * Init header offsets.
     *
     * We store in order:
     *     Header
     *     Samples
     *     Bases
     *     Comments
     *     Private
     */
    scf->header.samples_offset = (uint_4)sizeof(Header);
    size = scf->header.samples * (scf->header.sample_size == 1 ?
				  sizeof(Samples1) : sizeof(Samples2));
    scf->header.bases_offset = (uint_4)(scf->header.samples_offset +
					 size);
    size = scf->header.bases * sizeof(Bases);
    scf->header.comments_offset = (uint_4)(scf->header.bases_offset + size);

    size = scf->header.comments_size;
    scf->header.private_offset = (uint_4)(scf->header.comments_offset + size);

    /* Init a few other things, such as the magic number */
    scf->header.magic_number = SCF_MAGIC;

    if (scf_version == 3) {
	memcpy(scf->header.version, scf_version_float2str(SCF_VERSION), 4);
    } else {
	memcpy(scf->header.version, scf_version_float2str(SCF_VERSION_OLD), 4);
    }

    /* Write header */
    if (write_scf_header(fp, &scf->header) == -1)
	return -1;

    if (scf_version == 3) {
	/* Write Samples */
	if (scf->header.sample_size == 1)
	    err = write_scf_samples31(fp, scf->samples.samples1,
				      scf->header.samples);
	else
	    err = write_scf_samples32(fp, scf->samples.samples2,
				      scf->header.samples);
	if (-1 == err)
	    return -1;
	
	/* Write Bases */
	if (-1 == write_scf_bases3(fp, scf->bases, scf->header.bases))
	    return -1;

    } else {
	/* Write Samples */
	if (scf->header.sample_size == 1)
	    err = write_scf_samples1(fp, scf->samples.samples1,
				     scf->header.samples);
	else
	    err = write_scf_samples2(fp, scf->samples.samples2,
				     scf->header.samples);
	if (-1 == err)
	    return -1;
	
	/* Write Bases */
	if (-1 == write_scf_bases(fp, scf->bases, scf->header.bases))
	    return -1;
    }

    /* Write Comments */
    if (-1 == write_scf_comment(fp, scf->comments,
				scf->header.comments_size))
	return -1;

    /* Write private data */
    if (scf->header.private_size) {
	if (scf->header.private_size  != fwrite(scf->private_data, 1,
						scf->header.private_size, fp))
	    return -1;
    }

    return 0;
}

/*
 * Write Seq out as a .scf file to file 'fn'.
 */
int write_scf(Scf *scf, char *fn)
{
    FILE *fp;

    /* Open for for write in binary mode */
    if ((fp = fopen(fn,"wb")) == NULL) 
	return -1;

    if (fwrite_scf(scf, fp)) {
	fclose(fp);
	return -1;
    }

    fclose(fp);
    return 0;
}
