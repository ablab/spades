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
  Title:       read_scf.c
  
  Purpose:	 read IO of Standard Chromatogram Format sequences
  Last update:   August 18 1994
  
  Change log:
  4 Feb 1992,  Now draft proposal version 2
  20 Feb 1992, Grab info from comment lines
  19 Aug 1992, If SCF file has clip information, don't clip automatically
  10 Nov 1992  SCF comments now stored in seq data structure
  18 Aug 1994  Renamed from  ReadIOSCF.c; now purely SCF IO (no Seq structs)

*/

/* ---- Imports ---- */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <ctype.h>
#include <stdio.h>    /* IMPORT: fopen, fclose, fseek, ftell, fgetc,
			 EOF */
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "io_lib/mach-io.h"
#include "io_lib/xalloc.h"
#include "io_lib/compress.h"
#include "io_lib/Read.h"

#include "io_lib/stdio_hack.h"
#include "io_lib/scf.h"      /* SCF structures */


/* SunOS4 has it's definitions in unistd, which we won't include for compat. */
#ifndef SEEK_SET
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif

/* ---- Exported functions ---- */

int read_scf_header(FILE *fp, Header *h)
{
    int i;

    if (be_read_int_4(fp,&h->magic_number)==0)            return -1;

    if (h->magic_number != SCF_MAGIC)                     return -1;

    if (be_read_int_4(fp,&h->samples)==0)                 return -1;
    if (be_read_int_4(fp,&h->samples_offset)==0)          return -1;
    if (be_read_int_4(fp,&h->bases)==0)                   return -1;
    if (be_read_int_4(fp,&h->bases_left_clip)==0)         return -1;
    if (be_read_int_4(fp,&h->bases_right_clip)==0)        return -1;
    if (be_read_int_4(fp,&h->bases_offset)==0)            return -1;
    if (be_read_int_4(fp,&h->comments_size)==0)           return -1;
    if (be_read_int_4(fp,&h->comments_offset)==0)         return -1;
    if (fread(&h->version[0],sizeof(h->version),1,fp)!=1) return -1;
    if (be_read_int_4(fp,&h->sample_size)==0)             return -1;
    if (be_read_int_4(fp,&h->code_set)==0)                return -1;
    if (be_read_int_4(fp,&h->private_size)==0)            return -1;
    if (be_read_int_4(fp,&h->private_offset)==0)          return -1;
    for (i=0;i<18;i++)
	if (be_read_int_4(fp,&h->spare[i])==0)            return -1;
    
    return 0;
}


int read_scf_sample1(FILE *fp, Samples1 *s)
{
    uint_1 buf[4];

    if (4 != fread(buf, 1, 4, fp)) return -1;
    s->sample_A = buf[0];
    s->sample_C = buf[1];
    s->sample_G = buf[2];
    s->sample_T = buf[3];

/*
    if (1 != fread(s, 4, 1, fp)) return -1;
*/

    return 0;
}


int read_scf_sample2(FILE *fp, Samples2 *s)
{
    uint_2 buf[4];

    if (4 != fread(buf, 2, 4, fp)) return -1;
    s->sample_A = be_int2(buf[0]);
    s->sample_C = be_int2(buf[1]);
    s->sample_G = be_int2(buf[2]);
    s->sample_T = be_int2(buf[3]);
    
    return 0;
}

int read_scf_samples1(FILE *fp, Samples1 *s, size_t num_samples) {
    size_t i;

    for (i = 0; i < num_samples; i++) {
	if (-1 == read_scf_sample1(fp, &(s[i])))
	    return -1;
    }

    return 0;
}


int read_scf_samples2(FILE *fp, Samples2 *s, size_t num_samples) {
    size_t i;

    for (i = 0; i < num_samples; i++) {
	if (-1 == read_scf_sample2(fp, &(s[i])))
	    return -1;
    }

    return 0;
}


int read_scf_samples32(FILE *fp, Samples2 *s, size_t num_samples) {
    size_t i;
    uint2 *samples_out;

    /* version to read delta delta data in 2 bytes */

    if ( ! (samples_out = (uint2 *)xmalloc((num_samples+1) * 
					    sizeof(uint2)))) {
	return -1;
    }


    if (num_samples != fread(samples_out, 2, num_samples, fp)) return -1;
#ifdef SP_LITTLE_ENDIAN
    for (i = 0; i < num_samples; i++) {
	samples_out[i] = be_int2(samples_out[i]);
    }
#endif
    scf_delta_samples2 ( samples_out, num_samples, 0);
    for (i = 0; i < num_samples; i++) {
	(&s[i])->sample_A = samples_out[i];
    }

    if (num_samples != fread(samples_out, 2, num_samples, fp)) return -1;
#ifdef SP_LITTLE_ENDIAN
    for (i = 0; i < num_samples; i++) {
	samples_out[i] = be_int2(samples_out[i]);
    }
#endif
    scf_delta_samples2 ( samples_out, num_samples, 0);
    for (i = 0; i < num_samples; i++) {
	(&s[i])->sample_C = samples_out[i];
    }

    if (num_samples != fread(samples_out, 2, num_samples, fp)) return -1;
#ifdef SP_LITTLE_ENDIAN
    for (i = 0; i < num_samples; i++) {
	samples_out[i] = be_int2(samples_out[i]);
    }
#endif
    scf_delta_samples2 ( samples_out, num_samples, 0);
    for (i = 0; i < num_samples; i++) {
	(&s[i])->sample_G = samples_out[i];
    }

    if (num_samples != fread(samples_out, 2, num_samples, fp)) return -1;
#ifdef SP_LITTLE_ENDIAN
    for (i = 0; i < num_samples; i++) {
	samples_out[i] = be_int2(samples_out[i]);
    }
#endif
    scf_delta_samples2 ( samples_out, num_samples, 0);
    for (i = 0; i < num_samples; i++) {
	(&s[i])->sample_T = samples_out[i];
    }
    xfree(samples_out);
    return 0;
}

int read_scf_samples31(FILE *fp, Samples1 *s, size_t num_samples) {
    size_t i;
    int1 *samples_out;

    /* version to read delta delta data in 1 byte */

    if ( ! (samples_out = (int1 *)xmalloc((num_samples+1) * 
					    sizeof(int1)))) {
	return -1;
    }

    if (num_samples != fread(samples_out, 1, num_samples, fp)) return -1;
    scf_delta_samples1 ( samples_out, num_samples, 0);
    for (i = 0; i < num_samples; i++) {
	(&s[i])->sample_A = samples_out[i];
    }

    if (num_samples != fread(samples_out, 1, num_samples, fp)) return -1;
    scf_delta_samples1 ( samples_out, num_samples, 0);
    for (i = 0; i < num_samples; i++) {
	(&s[i])->sample_C = samples_out[i];
    }

    if (num_samples != fread(samples_out, 1, num_samples, fp)) return -1;
    scf_delta_samples1 ( samples_out, num_samples, 0);
    for (i = 0; i < num_samples; i++) {
	(&s[i])->sample_G = samples_out[i];
    }

    if (num_samples != fread(samples_out, 1, num_samples, fp)) return -1;
    scf_delta_samples1 ( samples_out, num_samples, 0);
    for (i = 0; i < num_samples; i++) {
	(&s[i])->sample_T = samples_out[i];
    }

    xfree(samples_out);
    return 0;
}

int read_scf_base(FILE *fp, Bases *b)
{
    uint_1 buf[12];

    if (1 != fread(buf, 12, 1, fp)) return -1;
    b->peak_index = be_int4(((uint_4 *)buf)[0]);
    b->prob_A = buf[4];
    b->prob_C = buf[5];
    b->prob_G = buf[6];
    b->prob_T = buf[7];
    b->base   = buf[8];
    b->spare[0] = buf[9];
    b->spare[1] = buf[10];
    b->spare[2] = buf[11];

    return 0;
}


int read_scf_bases(FILE *fp, Bases *b, size_t num_bases) {
    size_t i;

    for (i = 0; i < num_bases; i++) {
	if (-1 == read_scf_base(fp, &(b[i])))
	    return -1;
    }

    return 0;
}

int read_scf_bases3(FILE *fp, Bases *b, size_t num_bases)
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

    if (num_bases != fread(buf4, 4, num_bases, fp)) return -1;
    for (i=0; i < num_bases; i++)
	(&b[i])->peak_index = be_int4(buf4[i]);

    if (8 * num_bases != fread(buf1, 1, 8 * num_bases, fp)) return -1;

    for (i=0; i < num_bases; i++) {
	(&b[i])->prob_A   = buf1[i];
	(&b[i])->prob_C   = buf1[i+num_bases];
	(&b[i])->prob_G   = buf1[i+2*num_bases];
	(&b[i])->prob_T   = buf1[i+3*num_bases];
	(&b[i])->base     = buf1[i+4*num_bases];
	(&b[i])->spare[0] = buf1[i+5*num_bases];
	(&b[i])->spare[1] = buf1[i+6*num_bases];
	(&b[i])->spare[2] = buf1[i+7*num_bases];
    }

    xfree(buf4);
    xfree(buf1);

    return 0;
}



int read_scf_comment(FILE *fp, Comments *c, size_t s)
{
    if (fread(c, 1, s, fp) != s) return -1;

    return 0;
}


/*
 * Read the SCF format sequence from FILE *fp into a 'scf' structure.
 * A NULL result indicates failure.
 */
Scf *fread_scf(FILE *fp) {
    Scf *scf;
    Header h;
    int err;
    float scf_version;
    int sections = read_sections(0);

    /* Read header */
    if (read_scf_header(fp, &h) == -1) {
	return NULL;
    }

    /* Allocate memory */
    if (NULL == (scf = scf_allocate(h.samples, h.sample_size,
				    h.bases, h.comments_size,
				    h.private_size))) 
	return NULL;

    /* fake things for older style SCF -- SD */
    if (h.sample_size != 1 && h.sample_size != 2) h.sample_size = 1;

    scf_version = scf_version_str2float(h.version);

    memcpy(&scf->header, &h, sizeof(Header));

    if (sections & READ_SAMPLES) {
	/* Read samples */
	if (fseek(fp, (off_t)h.samples_offset, 0 /* SEEK_SET */) != 0) {
	    scf_deallocate(scf);
	    return NULL;
	}

	if ( 2.9 > scf_version ) {

	    if (h.sample_size == 1) {
		err= read_scf_samples1(fp, scf->samples.samples1, h.samples);
	    }
	    else {
		err= read_scf_samples2(fp, scf->samples.samples2, h.samples);
	    }
	}
	else {

	    if (h.sample_size == 1) {
		err= read_scf_samples31(fp, scf->samples.samples1, h.samples);
	    } 
	    else {
		err= read_scf_samples32(fp, scf->samples.samples2, h.samples);
	    }
	}
	if (-1 == err) {
	    scf_deallocate(scf);
	    return NULL;
	}
    }

    if (sections & READ_BASES) {
	/* Read bases */
	if (fseek(fp, (off_t)h.bases_offset, 0 /* SEEK_SET */) != 0) {
	    scf_deallocate(scf);
	    return NULL;
	}

	if ( 2.9 > scf_version ) {

	    if (-1 == read_scf_bases(fp, scf->bases, h.bases)) {
		scf_deallocate(scf);
		return NULL;
	    }
	}
	else {
	
	    if (-1 == read_scf_bases3(fp, scf->bases, h.bases)) {
		scf_deallocate(scf);
		return NULL;
	    }
	}
    }

    if (sections & READ_COMMENTS) {
	/* Read comments */
	if (scf->comments) {
	    if (fseek(fp,(off_t)(h.comments_offset), 0) != 0
		|| -1 == read_scf_comment(fp, scf->comments,
					  h.comments_size)) {
		/*
		 * Was: "scf_deallocate(scf); return NULL;".
		 * We now simply clear the comments and gracefully continue.
		 */
		fprintf(stderr, "Warning: SCF file had invalid comment field\n");
		xfree(scf->comments);
		scf->comments = NULL;
	    } else {
		scf->comments[h.comments_size] = '\0';
	    }
	}
    }

    /* Read private data */
    if (h.private_size) {
	if (-1 == fseek(fp, (off_t)(h.private_offset), 0) ||
	    h.private_size != fread(scf->private_data, 1, h.private_size, fp)){
	    scf_deallocate(scf);
	    return NULL;
	}
    }

    return scf;
}

/*
 * Read the SCF format sequence with name `fn' into a 'scf' structure.
 * A NULL result indicates failure.
 */
Scf *read_scf(char *fn) {
    Scf *scf;

    FILE *fp;

    /* Open fn for reading in binary mode */

    if (NULL == (fp = fopen_compressed(fn, NULL)))
	return NULL;

    scf = fread_scf(fp);
    fclose(fp);

    return scf;
}
