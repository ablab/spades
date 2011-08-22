/*
 * Copyright (c) Medical Research Council 1998. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies and that credit is given
 * where due.
 *
 * This file was written by James Bonfield, Simon Dear, Rodger Staden,
 * Kathryn Beal, as part of the Staden Package at the MRC Laboratory of
 * Molecular Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 */

/*
 * This file contains the necessary code for reading the quality values from
 * an SCF file. It supports both V2 and V3 SCF formats.
 * It's done in an efficient manner by extracting only the relevant SCF
 * components.
 * This file is derived from the Gap4 source file scf_extras.c.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdlib.h>

#include "io_lib/stdio_hack.h"
#include "io_lib/compress.h"
#include "io_lib/misc.h"
#include "io_lib/scf.h"
#include "io_lib/expFileIO.h"
#include "io_lib/traceType.h"
#include "io_lib/open_trace_file.h"
#include "io_lib/scf_extras.h"

/*
 * ---------------------------------------------------------------------------
 * Loads confidence values from the trace file and averages them.
 * 'opos' is optional - if not known then set to NULL.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int get_read_conf(Exp_info *e, int length, int2 *opos, int1 *conf) {
    int ttype, i;
    FILE *fp;
    uint_1 *prob_A, *prob_C, *prob_G, *prob_T;
    char *seq;
    float scf_version;
    int nbases = 0;

    /* Sanity check */
    if (!(exp_Nentries(e,EFLT_LT) && exp_Nentries(e,EFLT_LN)))
	return -1;

    /* Find and load trace file */
    ttype = trace_type_str2int(exp_get_entry(e, EFLT_LT));

    if (ttype != TT_SCF &&
	ttype != TT_CTF &&
	ttype != TT_ZTR)
	return -1;

    /*
     * We only support direct reading accuracy values from SCF files.
     * Otherwise we have to take a slower approach.
     */
    if (ttype != TT_SCF) {
	Read *r;
	int sec = read_sections(0);
	read_sections(READ_BASES);

	if (NULL == (r = read_reading(exp_get_entry(e,EFLT_LN), TT_ANYTR))) {
	    read_sections(sec);
	    return -1;
	}

	prob_A = (int1 *)xmalloc(r->NBases);
	prob_C = (int1 *)xmalloc(r->NBases);
	prob_G = (int1 *)xmalloc(r->NBases);
	prob_T = (int1 *)xmalloc(r->NBases);
	seq    = (char *)xmalloc(r->NBases);

	memcpy(prob_A, r->prob_A, r->NBases);
	memcpy(prob_C, r->prob_C, r->NBases);
	memcpy(prob_G, r->prob_G, r->NBases);
	memcpy(prob_T, r->prob_T, r->NBases);
	memcpy(seq,    r->base,   r->NBases);

	nbases = r->NBases;

	read_deallocate(r);
	read_sections(sec);

    } else {
	Header h;
	/* For SCF files we read directly - the above code would also do. */

	if (NULL == (fp = open_trace_file(exp_get_entry(e,EFLT_LN), NULL)))
	    return -1;

	/* Read the SCF header */
	if (-1 == read_scf_header(fp, &h))
	    return -1;
	scf_version = scf_version_str2float(h.version);
	nbases = h.bases;

	/* Alloc memory */
	prob_A = (uint_1 *)xmalloc(h.bases * sizeof(*prob_A));
	prob_C = (uint_1 *)xmalloc(h.bases * sizeof(*prob_A));
	prob_G = (uint_1 *)xmalloc(h.bases * sizeof(*prob_A));
	prob_T = (uint_1 *)xmalloc(h.bases * sizeof(*prob_A));
	seq    = (char   *)xmalloc(h.bases * sizeof(*seq));
	if (NULL == prob_A ||
	    NULL == prob_C ||
	    NULL == prob_G ||
	    NULL == prob_T ||
	    NULL == seq)
	    return -1;

	/* Load base scores */
	if (scf_version >= 3.0) {
	    /*
	     * Version 3 base format:
	     * num_bases * 4byte peak index
	     * num_bases * prob_A
	     * num_bases * prob_C
	     * num_bases * prob_G
	     * num_bases * prob_T
	     * num_bases * base
	     * num_bases * spare (x3)
	     */
	    fseek(fp, (off_t)h.bases_offset + 4 * h.bases, SEEK_SET);
	    if (h.bases != fread(prob_A, 1, h.bases, fp))
		return -1;
	    if (h.bases != fread(prob_C, 1, h.bases, fp))
		return -1;
	    if (h.bases != fread(prob_G, 1, h.bases, fp))
		return -1;
	    if (h.bases != fread(prob_T, 1, h.bases, fp))
		return -1;
	    if (h.bases != fread(seq, 1, h.bases, fp))
		return -1;
	} else {
	    int i;
	    uint_1 buf[12];

	    /*
	     * Version 2 base format
	     * num_bases * base_struct,  where base_struct is 12 bytes:
	     *     0-3 peak_index
	     *     4-7 prob_A/C/G/T
	     *     8   base
	     *     9-  spare
	     */
	    fseek(fp, (off_t)h.bases_offset, SEEK_SET);

	    for (i = 0; (unsigned)i < h.bases; i++) {
		if (1 != fread(buf, 12, 1, fp))
		    return -1;
		prob_A[i] = buf[4];
		prob_C[i] = buf[5];
		prob_G[i] = buf[6];
		prob_T[i] = buf[7];
		seq[i]    = buf[8];
	    }
	}

	fclose(fp);
    }

    /* Determine confidence values */
    if (opos) {
	for (i=0; i<length; i++) {
	    if (opos[i] == 0) {
		/* Inserted base, change to 0% */
		conf[i] = 0;
	    } else {
		switch(seq[opos[i]-1]) {
		case 'a':
		case 'A':
		    conf[i] = prob_A[opos[i]-1];
		    break;
		case 'c':
		case 'C':
		    conf[i] = prob_C[opos[i]-1];
		    break;
		case 'g':
		case 'G':
		    conf[i] = prob_G[opos[i]-1];
		    break;
		case 't':
		case 'T':
		    conf[i] = prob_T[opos[i]-1];
		    break;
		default:
		    conf[i] = 2;
		}
	    }
	}
    } else {
	int mlength = MIN(length, nbases);

	for (i=0; i < mlength; i++) {
	    switch(seq[i]) {
	    case 'a':
	    case 'A':
		conf[i] = prob_A[i];
		break;
	    case 'c':
	    case 'C':
		conf[i] = prob_C[i];
		break;
	    case 'g':
	    case 'G':
		conf[i] = prob_G[i];
		break;
	    case 't':
	    case 'T':
		conf[i] = prob_T[i];
		break;
	    case 'n':
	    case 'N':
	    case '-':
		conf[i] = (prob_A[i] + prob_C[i] + prob_G[i] + prob_T[i]) / 4;
		break;
	    default:
		conf[i] = 2;
	    }
	}
	for (; i < length; i++)
	    conf[i] = 2;
    }

    xfree(prob_A);
    xfree(prob_C);
    xfree(prob_G);
    xfree(prob_T);
    xfree(seq);

    return 0;
}
