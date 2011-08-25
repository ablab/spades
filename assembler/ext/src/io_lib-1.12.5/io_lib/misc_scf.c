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
  Title:       misc_scf.c
  
  Purpose:	 misc handling of Standard Chromatogram Format sequences
  Last update:   August 18 1994
  
  Change log:
  18 Aug 1994  Creation from bits of {read,write}scf.c and new code.

*/

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <string.h>

#include "io_lib/scf.h"
#include "io_lib/mach-io.h"
#include "io_lib/xalloc.h"

#include "io_lib/stdio_hack.h"

float scf_version_str2float(char version[])
{
    char v[5];
    strncpy(v,version,4);v[4]='\0';
    if (strspn(v,"0123456789. ")!=4) return 0.0;
    return (float)atof(v);
}

char *scf_version_float2str(float f)
{
    static char v[5];

    sprintf(v, "%1.2f", f);
    return v;
}


/*
 * Allocates memory for the scf elements based upon arguments passed.
 * Returns;
 *    Scf *	- Success. The scf structure and it's samples, bases,
 *                and comments fields have been allocated.
 *    NULL	- Failure.
 */
Scf *scf_allocate(int num_samples, int sample_size, int num_bases,
		  int comment_size, int private_size) {
    Scf *scf;

    scf = (Scf *)xcalloc(1, sizeof(Scf));
    if (NULL == scf)
	return NULL;

    /* bases - +1 as a safety guard for when num_bases==0 */
    scf->bases = (Bases *)xcalloc(sizeof(Bases), num_bases+1);
    if (NULL == scf->bases)
	return NULL;

    /* samples */
    scf->header.sample_size = sample_size;
    if (scf->header.sample_size == 1) {
	scf->samples.samples1 = (Samples1 *)xmalloc(num_samples *
						    sizeof(Samples1) + 1);
    } else {
	scf->samples.samples2 = (Samples2 *)xmalloc(num_samples *
						    sizeof(Samples2) + 1);
    }
    if (NULL == scf->samples.samples1) {
	xfree(scf->bases);
	xfree(scf);
	return NULL;
    }

    /* comments */
    if (comment_size) {
	scf->comments = (Comments *)xmalloc(sizeof(Comments) *
					    (comment_size + 1));
	if (NULL == scf->comments) {
	    xfree(scf->bases);
	    xfree(scf->samples.samples1);
	    xfree(scf);
	    return NULL;
	}
    } else
	scf->comments = NULL;

    /* private data */
    if (private_size) {
	scf->private_data = (char *)xmalloc(private_size);
	if (NULL == scf->private_data) {
	    xfree(scf->bases);
	    xfree(scf->samples.samples1);
	    if (scf->comments) xfree(scf->comments);
	    xfree(scf);
	    return NULL;
	}
    } else
	scf->private_data = NULL;
    
    return scf;
}

void scf_deallocate(Scf *scf) {
    xfree(scf->bases);
    xfree(scf->samples.samples1);
    if (scf->comments)
	xfree(scf->comments);
    if (scf->private_data)
	xfree(scf->private_data);
    xfree(scf);
}


int is_scf(char *fn)
/*
 * Check to see if file with name `fn' is in SCF format
 * 
 * Returns:
 * 1 - is SCF format
 * 0 - not SCF format
 *-1 - failure
 */
{
    FILE *fp;
    uint_4 magic;
    int ok;
    
    if ( (fp=fopen(fn,"rb")) == NULL) {
	ok = -1;
    } else {
	if ( be_read_int_4(fp, &magic) != 1 ) {
	    ok = 0;
	} else {
	    ok = (magic==SCF_MAGIC);
	}
	fclose(fp);
    }
    
    return ok;
}

void scf_delta_samples1 ( int1 samples[], int num_samples, int job) {

    /* If job == DELTA_IT:
       change a series of sample points to a series of delta delta values:
       ie change them first: delta = current_value - previous_value
       then delta_delta = delta - previous_delta

       else
       do the reverse
       */

    int i;

    if ( DELTA_IT == job ) {
#ifdef CLEAR_BUT_SLOW
	int1 p_delta, p_sample;

	p_delta  = 0;
	for (i=0;i<num_samples;i++) {
	    p_sample = samples[i];
	    samples[i] = samples[i] - p_delta;
	    p_delta  = p_sample;
	}
	p_delta  = 0;
	for (i=0;i<num_samples;i++) {
	    p_sample = samples[i];
	    samples[i] = samples[i] - p_delta;
	    p_delta  = p_sample;
	}
#else
	for (i = num_samples-1 ; i > 1; i--) {
	    samples[i] = samples[i] - 2*samples[i-1] + samples[i-2];
	}
	samples[1] = samples[1] - 2*samples[0];
#endif

    } else {

#ifdef CLEAR_BUT_SLOW
	int1 p_sample;

	p_sample = 0;
	for (i=0;i<num_samples;i++) {
	    samples[i] = samples[i] + p_sample;
	    p_sample = samples[i];
	}
	p_sample = 0;
	for (i=0;i<num_samples;i++) {
	    samples[i] = samples[i] + p_sample;
	    p_sample = samples[i];
	}
#else
	int1 p_sample1, p_sample2;
	
	p_sample1 = p_sample2 = 0;
	for (i = 0; i < num_samples; i++) {
	    p_sample1  = p_sample1 + samples[i];
	    samples[i] = p_sample1 + p_sample2;
	    p_sample2  = samples[i];
	}
#endif
    }
}

void scf_delta_samples2 ( uint_2 samples[], int num_samples, int job) {

    /* If job == DELTA_IT:
       change a series of sample points to a series of delta delta values:
       ie change them first: delta = current_value - previous_value
       then delta_delta = delta - previous_delta

       else
       do the reverse
       */

    register int i;

    if ( DELTA_IT == job ) {
#ifdef CLEAR_BUT_SLOW
	register uint_2 p_delta, p_sample;

	p_delta  = 0;
	for (i=0;i<num_samples;i++) {
	    p_sample = samples[i];
	    samples[i] = samples[i] - p_delta;
	    p_delta  = p_sample;
	}
	p_delta  = 0;
	for (i=0;i<num_samples;i++) {
	    p_sample = samples[i];
	    samples[i] = samples[i] - p_delta;
	    p_delta  = p_sample;
	}
#else
	for (i = num_samples-1 ; i > 1; i--) {
	    samples[i] = samples[i] - 2*samples[i-1] + samples[i-2];
	}
	samples[1] = samples[1] - 2*samples[0];
#endif

    } else {

#ifdef CLEAR_BUT_SLOW
	register uint_2 p_sample;

	p_sample = 0;
	for (i=0;i<num_samples;i++) {
	    samples[i] = samples[i] + p_sample;
	    p_sample = samples[i];
	}
	p_sample = 0;
	for (i=0;i<num_samples;i++) {
	    samples[i] = samples[i] + p_sample;
	    p_sample = samples[i];
	}
#else
	uint_2 p_sample1, p_sample2;
	
	p_sample1 = p_sample2 = 0;
	for (i = 0; i < num_samples; i++) {
	    p_sample1  = p_sample1 + samples[i];
	    samples[i] = p_sample1 + p_sample2;
	    p_sample2  = samples[i];
	}
#endif
    }
}

