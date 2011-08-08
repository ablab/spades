/*
 * Copyright (c) Medical Research Council 1994-1998. All rights reserved.
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
 * makeSCF v3.06, 17/04/2001
 *
 * Derived from the older makeSCF; this one has been rewritten to use the new
 * IO libraries and no longer performs quality clipping itself. It also writes
 * in a new format that is more easily compressed.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <strings.h>
#include <io_lib/Read.h>
#include <io_lib/traceType.h>
#include <io_lib/xalloc.h>

/*
 * Add our comments.
 * 1. Work out the comment length - simply add 1K to the current length.
 * 2. Copy the old comments to our new.
 * 3. Replace '\n' with ' '.
 * 4. Add our comments and switch.
 */
void add_comments(Read *r, char *name, int format) {
    Comments *cc;
    int clen;
    char *cp;
    char buf[1024];

    /* 1. */
    if (r->info) {
	clen = strlen(r->info) + 1024;
    } else {
	clen = 1024;
    }

    /* 2. */
    if (NULL == (cc = (char *)xmalloc(clen)))
	return;

    if (r->info)
	strcpy(cc, r->info);
    else
	*cc = 0;
    
    /* 3. */
    cp = cc;
/*
    while (*cp) {
	if (*cp == '\n')
	    *cp = ' ';
	
	cp++;
    }
    *cp++ = '\n';
*/
  
    /* 4. */
    sprintf(buf, "CONV=makeSCF V3.06\nDATF=%s\nDATN=%s\n",
	    trace_type_int2str(format), name);

    strcat(cp, buf);

    if (r->info)
	xfree(r->info);

    r->info = cc;
}

void scale_trace8(Read *r) {
    double s;
    int i;

    if (r->maxTraceVal <= 255)
	return;
    
    s = ((double)255)/r->maxTraceVal;

    for (i = 0; i < r->NPoints; i++) {
	r->traceA[i] *= s;
	r->traceC[i] *= s;
	r->traceG[i] *= s;
	r->traceT[i] *= s;
    }

    r->maxTraceVal = 255;
}

/*
 * Here we just take the minimum trace value and subtract this from all others.
 * The assumption is that the signal will always be 'base line' on at least
 * one of the four channels.
 */
void subtract_background(Read *r) {
    int i, min;
    for (i = 0; i < r->NPoints; i++) {
	min = 999999;
	if (r->traceA[i] < min) min = r->traceA[i];
	if (r->traceC[i] < min) min = r->traceC[i];
	if (r->traceG[i] < min) min = r->traceG[i];
	if (r->traceT[i] < min) min = r->traceT[i];
	r->traceA[i] -= min;
	r->traceC[i] -= min;
	r->traceG[i] -= min;
	r->traceT[i] -= min;
    }
}

/*
 * Find the average background level of a trace, and subtract this from the
 * peak heights.
 *
 * NB. That method is flawed. For now take the minimum instead of average, but
 * this also has horrid flaws. See above method.
 */
void subtract_background_old(Read *r) {
    int i, j, min, bg, max = 0;
    int win_len = 501, win_len2 = win_len/2;
    int *background;

    if (NULL == (background = (int *)xmalloc((r->NPoints + 2 * win_len)
					     * sizeof(*background))))
	return;

    if (r->NPoints < win_len)
	win_len = r->NPoints;

    /* Find minimum trace levels at each point */
    for (i = 0, j = win_len2; i < r->NPoints; i++, j++) {
	min = INT_MAX;
	if (r->traceA[i] < min)
	    min = r->traceA[i];
	if (r->traceC[i] < min)
	    min = r->traceC[i];
	if (r->traceG[i] < min)
	    min = r->traceG[i];
	if (r->traceT[i] < min)
	    min = r->traceT[i];
	background[j] = min;
    }
    for (i = 0; i < win_len2; i++) {
	background[i] = background[i + win_len2];
	background[i + r->NPoints] = background[i + r->NPoints - win_len2];
    }

    /* Take lowest background over win_len and subtract it */
    for (i = 0; i < r->NPoints; i++) {
	/* Could optimise this considerably */
	bg = INT_MAX;
	for (j = 0; j < win_len; j++) {
	    if (background[i + j] < bg)
		bg = background[i + j];
	}

	r->traceA[i] -= bg;
	r->traceC[i] -= bg;
	r->traceG[i] -= bg;
	r->traceT[i] -= bg;

	if (r->traceA[i] > max) max = r->traceA[i];
	if (r->traceC[i] > max) max = r->traceC[i];
	if (r->traceG[i] > max) max = r->traceG[i];
	if (r->traceT[i] > max) max = r->traceT[i];
    }
    
    r->maxTraceVal = max;

    xfree(background);
}

/*
 * Find the maximum height of traces at the called bases. Use this to clip any
 * other bases.
 */
void reset_max_called_height(Read *r) {
    int i, max = 0;

    /* Find max */
    for (i=0; i < r->NBases; i++) {
	switch(r->base[i]) {
	case 'a':
	case 'A':
	    if (r->traceA[r->basePos[i]] > max)
		max = r->traceA[r->basePos[i]];
	    break;

	case 'c':
	case 'C':
	    if (r->traceC[r->basePos[i]] > max)
		max = r->traceC[r->basePos[i]];
	    break;

	case 'g':
	case 'G':
	    if (r->traceG[r->basePos[i]] > max)
		max = r->traceG[r->basePos[i]];
	    break;

	case 't':
	case 'T':
	    if (r->traceT[r->basePos[i]] > max)
		max = r->traceT[r->basePos[i]];
	    break;
	}
    }

    /* Clip to max */
    for (i = 0; i < r->NPoints; i++) {
	if (r->traceA[i] > max)
	    r->traceA[i] = max;
	if (r->traceC[i] > max)
	    r->traceC[i] = max;
	if (r->traceG[i] > max)
	    r->traceG[i] = max;
	if (r->traceT[i] > max)
	    r->traceT[i] = max;
    }
    if (r->maxTraceVal > max)
	r->maxTraceVal = max;
}

void rescale_heights(Read *r) {
    int win_len = 1000;
    int total = 0;
    int max, max2;
    int i, j, k;
    double max_val = 0, rescale = 1.0;
    TRACE *ta, *tc, *tg, *tt;

    ta = r->traceA;
    tc = r->traceC;
    tg = r->traceG;
    tt = r->traceT;

    if (r->NPoints < 2*win_len + 1)
	return;

    for (k = 0; k < 2; k++) {
	max2 = win_len * r->maxTraceVal;
	
	for (i = 0; i < win_len; i++) {
	    max = 0;
	    if (ta[i] > max) max = ta[i];
	    if (tc[i] > max) max = tc[i];
	    if (tg[i] > max) max = tg[i];
	    if (tt[i] > max) max = tt[i];
	    total += max;
	}
	
	for (j = 0; i < r->NPoints; i++, j++) {
	    max = 0;
	    if (ta[j] > max) max = ta[j];
	    if (tc[j] > max) max = tc[j];
	    if (tg[j] > max) max = tg[j];
	    if (tt[j] > max) max = tt[j];
	    total -= max;
	    
	    max = 0;
	    if (ta[i] > max) max = ta[i];
	    if (tc[i] > max) max = tc[i];
	    if (tg[i] > max) max = tg[i];
	    if (tt[i] > max) max = tt[i];
	    total += max;
	    
	    if (k == 0) {
		if (r->traceA[j] * ((double)max2 / total) > max_val)
		    max_val = r->traceA[j] * ((double)max2 / total);
		if (r->traceC[j] * ((double)max2 / total) > max_val)
		    max_val = r->traceC[j] * ((double)max2 / total);
		if (r->traceG[j] * ((double)max2 / total) > max_val)
		    max_val = r->traceG[j] * ((double)max2 / total);
		if (r->traceT[j] * ((double)max2 / total) > max_val)
		    max_val = r->traceT[j] * ((double)max2 / total);
	    } else {
		r->traceA[j] *= (double)max2 / total * rescale;
		r->traceC[j] *= (double)max2 / total * rescale;
		r->traceG[j] *= (double)max2 / total * rescale;
		r->traceT[j] *= (double)max2 / total * rescale;
	    }
	}

	for (; j < r->NPoints; j++) {
	    if (k == 0) {
		if (r->traceA[j] * ((double)max2 / total) > max_val)
		    max_val = r->traceA[j] * ((double)max2 / total);
		if (r->traceC[j] * ((double)max2 / total) > max_val)
		    max_val = r->traceC[j] * ((double)max2 / total);
		if (r->traceG[j] * ((double)max2 / total) > max_val)
		    max_val = r->traceG[j] * ((double)max2 / total);
		if (r->traceT[j] * ((double)max2 / total) > max_val)
		    max_val = r->traceT[j] * ((double)max2 / total);
	    } else {
		r->traceA[j] *= (double)max2 / total;
		r->traceC[j] *= (double)max2 / total;
		r->traceG[j] *= (double)max2 / total;
		r->traceT[j] *= (double)max2 / total;
	    }
	}

	if (max_val > 65535)
	    rescale = 65535 / max_val;
	else
	    rescale = 1.0;
    }
}

static int convert(char *in, mFILE *ofp, char *out, int format, int prec,
		   int comp, int normalise) {
    Read *r;

    if (NULL == (r = read_reading(in, format))) {
	fprintf(stderr, "%s: failed to read\n", in);
	return 1;
    }

    if (normalise) {
	subtract_background(r);
	reset_max_called_height(r);
	rescale_heights(r);
    }

    add_comments(r, in, format);
    if (prec == 1)
	scale_trace8(r);

    if (comp != -1)
	set_compression_method(comp);
    if (0 != (mfwrite_reading(ofp, r, TT_SCF))) {
	fprintf(stderr, "%s: failed to write\n", out);
	read_deallocate(r);
	return 1;
    }

    read_deallocate(r);
    return 0;
}


void usage(void) {
    fprintf(stderr,
	    "makeSCF [-8] [-2] [-3] [-s] [-compress mode] [-normalise]\n"
	    "       -(abi|alf|scf|any) input_name [-output output_name]\n"
	    " or\n"
	    "makeSCF [-8] [-2] [-3] [-s] [-compress mode] [-normalise]\n"
	    "       [-(abi|alf|scf|any)] input_name1 output_name1 ... "
	    "input_nameN output_nameN \n");

    exit(1);
}

int main(int argc, char **argv) {
    int format = TT_ANY, r, prec = 0, version = 3, silent = 0;
    int compress_mode = -1;
    char *inf = NULL;
    char *outf = NULL;
    mFILE *ofp = mstdout();
    int normalise = 0;

    for (argc--, argv++; argc > 0; argc--, argv++) {
	if (strcmp(*argv, "-8") == 0) {
	    prec = 1;
	} else if (strcmp(*argv, "-2") == 0) {
	    version = 2;
	} else if (strcmp(*argv, "-3") == 0) {
	    version = 3;
	} else if (strcmp(*argv, "-normalise") == 0) {
	    normalise = 1;
	} else if (strcmp(*argv, "-s") == 0) {
	    silent = 1;
	} else if (strcasecmp(*argv, "-abi") == 0) {
	    format = TT_ABI;
	    inf = *++argv;
	    argc--;
	} else if (strcasecmp(*argv, "-alf") == 0) {
	    format = TT_ALF;
	    inf = *++argv;
	    argc--;
	} else if (strcasecmp(*argv, "-scf") == 0) {
	    format = TT_SCF;
	    inf = *++argv;
	    argc--;
	} else if (strcasecmp(*argv, "-ztr") == 0) {
	    format = TT_ZTR;
	    inf = *++argv;
	    argc--;
	} else if (strcasecmp(*argv, "-any") == 0) {
	    format = TT_ANY;
	    inf = *++argv;
	    argc--;
	} else if (strcasecmp(*argv, "-output") == 0) {
	    outf = *++argv;
	    argc--;
	} else if (strcasecmp(*argv, "-compress") == 0) {
	    compress_mode = compress_str2int(*++argv);
	    argc--;
	} else {
	    break;
	}
    }

    /* if no args left than input file must have been specified */
    if (!argc && !inf)
	usage();

    /* if outfile set, then using original syntax, so don't expect
       any extra args */
    if (argc && outf)
      usage();

    if (!silent) {
	printf("makeSCF v3.06\n");
	printf("Copyright (c) MRC Laboratory of Molecular Biology, 2001. All rights reserved.\n");
    }

    set_scf_version(version);


    if(!argc) {
	/* original calling syntax */
	if (outf) {
	    ofp = mfopen(outf, "wb+");
	    if (NULL == ofp) {
		perror(outf);
		return 1;
	    }
	}

	r = convert(inf, ofp, outf, format, prec, compress_mode, normalise);
	mfclose(ofp);

	return r;

    }

    /* else */ {
	/* new calling syntax, handling multiple files */
	int result=0;

	for (; argc > 0; argc--, argv++) {
	    if (inf) {
		/* got infile, so get outfile and process */
		outf= *argv;
		ofp = mfopen(outf, "wb+");
		
		if (NULL == ofp) {
		    perror(outf);
		    if(!result) result=1;
		    continue;
		}
		r = convert(inf, ofp, outf, format, prec, compress_mode,
			    normalise);
		mfclose(ofp);
		if(!result) /* keep track of the first error */
		    result=r;
	      
		/* now need to get another infile */
		inf=NULL;
	    } else {
		/* need infile */
		inf= *argv;
	    }
	}

	return result;
    }
}

