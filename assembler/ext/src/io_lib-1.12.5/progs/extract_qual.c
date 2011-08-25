/*
 * Copyright (c) Medical Research Council 1994-1999. All rights reserved.
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

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <unistd.h>
#include <io_lib/Read.h>
#include <io_lib/traceType.h>
#include <io_lib/expFileIO.h>
#include <io_lib/open_trace_file.h>
#include <io_lib/xalloc.h>

/* #include "stdio_hack.h" */

#define LINE_LENGTH 60

/*
 * Converts the confidence array to the accuracy value string (AV).
 *
 * Note no memory overrun checks are performed on buf. It is recommended
 * that it is allocated to 4*len (worst case of "100 " for each base).
 *
 * Returns the buf argument.
 */
static char *my_conf2str(int1 *conf, int len, char *buf) {
    int i;
    char *ret = buf, *rs = buf;

    for (i = 0; i < len; i++) {
	sprintf(buf, "%d ", conf[i]);
	buf += strlen(buf);

	if (buf - rs > LINE_LENGTH) {
	    *buf++ = '\n';
	    *buf = '\0';
	    rs = buf;
	}
    }

    return ret;
}

static int do_trans(mFILE *infp, char *in_file, FILE *outfp, int format,
		    int good_only, int clip_cosmid, int fasta_out) {
    Read *r;
    char *tmp_prob_A, *tmp_prob_C, *tmp_prob_G, *tmp_prob_T;
    char *cstr = NULL;

    read_sections(READ_BASES);
    if (NULL == (r = mfread_reading(infp, in_file, format))) {
	fprintf(stderr, "Failed to read file '%s'\n", in_file);
	return 1;
    }

    tmp_prob_A = r->prob_A;
    tmp_prob_C = r->prob_C;
    tmp_prob_G = r->prob_G;
    tmp_prob_T = r->prob_T;

#ifdef IOLIB_EXP
    if (good_only && r->orig_trace_format == TT_EXP) {
	int left=0, right=r->NBases + 1, val, lval, rval;
	Exp_info *e = (Exp_info *)r->orig_trace;

	if (0 == exp_get_int(e, EFLT_SL, &val))
	    if (val > left)
		left = val;
	if (0 == exp_get_int(e, EFLT_QL, &val))
	    if (val > left)
		left = val;

	if (0 == exp_get_int(e, EFLT_SR, &val))
	    if (val < right)
		right = val;
	if (0 == exp_get_int(e, EFLT_QR, &val))
	    if (val < right)
		right = val;

	/* This is horrid - see gap seqInfo.c file for explaination */
	if (clip_cosmid) {
	    int got_cosmid;

	    if (0 == exp_get_rng(e, EFLT_CS, &lval, &rval)) {
		got_cosmid = 1;
	    } else if (0 == exp_get_int(e, EFLT_CL, &lval) &&
		       0 == exp_get_int(e, EFLT_CR, &rval)) {
		got_cosmid = 1;
	    } else {
		got_cosmid = 0;
	    }
	
	    if (got_cosmid) {
		if      (lval <= left   && rval <= left)    ;
		else if (lval <= left+1 && rval <  right)   left  = rval;
		else if (lval <= left+1 && rval >= right)   right = left+1;
		else if (lval <  right  && rval <  right)   right = lval;
		else if (lval <  right  && rval >= right)   right = lval;
	    }
	}

        r->prob_A += left;
        r->prob_C += left;
        r->prob_G += left;
        r->prob_T += left;
	r->NBases = right - left - 1;
    } else
#endif /* IOLIB_EXP */
    if (good_only) {
        r->prob_A += r->leftCutoff;
        r->prob_C += r->leftCutoff;
        r->prob_G += r->leftCutoff;
        r->prob_T += r->leftCutoff;
	r->NBases = r->rightCutoff - r->leftCutoff - 1;
    }

    /* Confidence values */
    if (r->prob_A && r->prob_C && r->prob_G && r->prob_T &&
	r->NBases > 0) {
	int i;

	/* We have some, but are they non zero values? */
        for (i = 0; i < r->NBases; i++) {
	    if (r->prob_A[i] || r->prob_C[i] ||
		r->prob_G[i] || r->prob_T[i])
		break;
	}
	if (i != r->NBases) {
	    int1 *conf = (int1 *)xmalloc(r->NBases);
	    cstr = (char *)xmalloc(r->NBases * 4 + 2);

	    for (i = 0; i < r->NBases; i++) {
		switch (r->base[i]) {
		case 'a':
		case 'A':
		    conf[i] = r->prob_A[i];
		    break;
		case 'c':
		case 'C':
		    conf[i] = r->prob_C[i];
		    break;
		case 'g':
		case 'G':
		    conf[i] = r->prob_G[i];
		    break;
		case 't':
		case 'T':
		    conf[i] = r->prob_T[i];
		    break;
                case 'b':
                case 'B':
		    conf[i] = (r->prob_C[i] + r->prob_G[i] + r->prob_T[i]) / 3;
                    break;
                case 'd':
                case 'D':
		    conf[i] = (r->prob_A[i] + r->prob_G[i] + r->prob_T[i]) / 3;
                    break;
                case 'h':
                case 'H':
		    conf[i] = (r->prob_A[i] + r->prob_C[i] + r->prob_T[i]) / 3;
                    break;
                case 'v':
                case 'V':
		    conf[i] = (r->prob_A[i] + r->prob_C[i] + r->prob_G[i]) / 3;
                    break;
                case 'k':
                case 'K':
		    conf[i] = (r->prob_G[i] + r->prob_T[i]) / 2;
                    break;
                case 'm':
                case 'M':
		    conf[i] = (r->prob_A[i] + r->prob_C[i]) / 2;
                    break;
                case 'r':
                case 'R':
		    conf[i] = (r->prob_A[i] + r->prob_G[i]) / 2;
                    break;
                case 's':
                case 'S':
		    conf[i] = (r->prob_C[i] + r->prob_G[i]) / 2;
                    break;
                case 'w':
                case 'W':
		    conf[i] = (r->prob_A[i] + r->prob_T[i]) / 2;
                    break;
                case 'y':
                case 'Y':
		    conf[i] = (r->prob_C[i] + r->prob_T[i]) / 2;
                    break;
		default:
		    conf[i] = (r->prob_A[i] + r->prob_C[i] + r->prob_G[i] + r->prob_T[i]) / 4;
		}
	    }

	    my_conf2str(conf, r->NBases, cstr);
	    xfree(conf);
	}
    }

    if (fasta_out) {
	char *p = strrchr(in_file, '/');
	/* Add header */
	if (NULL == p)
	    p = in_file;
	else
	    p++;
	fprintf(outfp, ">%s\n", p); 
    }

    if (cstr) {
      fprintf(outfp,"%s\n", cstr);
      xfree(cstr);
    }

    r->prob_A = tmp_prob_A;
    r->prob_C = tmp_prob_C;
    r->prob_G = tmp_prob_G;
    r->prob_T = tmp_prob_T;
    read_deallocate(r);
    fflush(outfp);

    return 0;
}

static void usage(void) {
    fprintf(stderr, "Usage: extract_qual [-r] [-(abi|alf|scf|exp|pln|ctf|ztr)]\n"
	    "                   [-good_only] [-clip_cosmid] [-fasta_out]\n"
	    "                   [-output output_name] [input_name] ...\n");
    exit(1);
}

int main(int argc, char **argv) {
    int from_stdin = 1;
    mFILE *infp = mstdin();
    FILE *outfp = stdout;
    int format = TT_ANY;
    int redirect = 1;
    int good_only = 0;
    int clip_cosmid = 0;
    int fasta_out = 0;
    int ret = 0;
    char *fofn = NULL;

    for (argc--, argv++; argc > 0; argc--, argv++) {
	if (strcmp(*argv, "-r") == 0) {
	    redirect = 0;
	} else if (strcasecmp(*argv, "-abi") == 0) {
            format = TT_ABI;
        } else if (strcasecmp(*argv, "-alf") == 0) {
            format = TT_ALF;
        } else if (strcasecmp(*argv, "-scf") == 0) {
            format = TT_SCF;
        } else if (strcasecmp(*argv, "-exp") == 0) {
            format = TT_EXP;
        } else if (strcasecmp(*argv, "-pln") == 0) {
            format = TT_PLN;
        } else if (strcasecmp(*argv, "-ztr") == 0) {
            format = TT_ZTR;
        } else if (strcasecmp(*argv, "-ctf") == 0) {
            format = TT_CTF;
        } else if (strcasecmp(*argv, "-good_only") == 0) {
	    good_only = 1;
        } else if (strcasecmp(*argv, "-clip_cosmid") == 0) {
	    clip_cosmid = 1;
	} else if (strcasecmp(*argv, "-fasta_out") == 0) {
	    fasta_out = 1;
	} else if (strcmp(*argv, "-fofn") == 0) {
	    fofn = *++argv;
	    argc--;
	    from_stdin = 0;
        } else if (strcasecmp(*argv, "-output") == 0) {
	    if (NULL == (outfp = fopen(*++argv, "wb"))) {
		perror(*argv);
		return 1;
	    }
            argc--;
	} else if (**argv != '-') {
	    from_stdin = 0;
	    break;
        } else {
            usage();
        }
    }

    read_experiment_redirect(redirect);

    if (!from_stdin) {
	if (fofn) {
	    FILE *fofn_fp;
	    char line[8192];

	    if (strcmp(fofn, "stdin") == 0)
		fofn_fp = stdin;
	    else
		fofn_fp = fopen(fofn, "r");

	    if (fofn_fp) {
		while (fgets(line, 8192, fofn_fp) != NULL) {
		    char *cp;
		    if (cp = strchr(line, '\n'))
			*cp = 0;
		    if (format == TT_EXP) {
			infp = open_exp_mfile(line, NULL);
		    } else {
			infp = open_trace_mfile(line, NULL);
		    }
		    if (NULL == infp) {
			perror(line);
			ret = 1;
		    } else {
			ret |= do_trans(infp, line, outfp, format, good_only,
					clip_cosmid, fasta_out);
			mfclose(infp);
		    }
		}
		fclose(fofn_fp);
	    }
	}
	for (;argc > 0; argc--, argv++) {
	    if (format == TT_EXP) {
		infp = open_exp_mfile(*argv, NULL);
	    } else {
		infp = open_trace_mfile(*argv, NULL);
	    }
	    if (NULL == infp) {
		perror(*argv);
		ret = 1;
	    } else {
		ret |= do_trans(infp, *argv, outfp, format, good_only,
				clip_cosmid, fasta_out);
		mfclose(infp);
	    }
	}
    } else {
	ret = do_trans(infp, "<stdin>", outfp, format, good_only, clip_cosmid,
		       fasta_out);
    }

    return ret;
}

