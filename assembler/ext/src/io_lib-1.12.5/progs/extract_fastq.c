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

static int do_trans(mFILE *infp, char *in_file, FILE *outfp, int format) {
    Read *r;
    char *p = strrchr(in_file, '/');
    int i;

    read_sections(READ_BASES);
    if (NULL == (r = mfread_reading(infp, in_file, format))) {
	fprintf(stderr, "Failed to read file '%s'\n", in_file);
	return 1;
    }

    if (NULL == p)
	p = in_file;
    else
	p++;

    fprintf(outfp, "@%s\n", p);
    fprintf(outfp, "%.*s\n", r->NBases, r->base);
    fprintf(outfp, "+%s\n", p);
    for (i = 0; i < r->NBases; i++) {
	int qual;
	switch (r->base[i]) {
	case 'A':
	case 'a':
	    qual = r->prob_A[i];
	    break;
	case 'C':
	case 'c':
	    qual = r->prob_C[i];
	    break;
	case 'G':
	case 'g':
	    qual = r->prob_G[i];
	    break;
	case 'T':
	case 't':
	    qual = r->prob_T[i];
	    break;
	default:
	    qual = 0;
	}
	fputc(qual + 33, outfp);
    }
    fputc('\n', outfp);

    read_deallocate(r);
    fflush(outfp);

    return 0;
}

static void usage(void) {
    fprintf(stderr, "Usage: extract_fastq [-(abi|alf|scf|exp|pln)]\n"
	    "                   [-output output_name] [-fofn fofn] [input_name] ...\n");
    exit(1);
}

int main(int argc, char **argv) {
    int from_stdin = 1;
    mFILE *infp = mstdin();
    FILE *outfp = stdout;
    int format = TT_ANY;
    int ret = 0;
    char *fofn = NULL;

    for (argc--, argv++; argc > 0; argc--, argv++) {
	if (strcasecmp(*argv, "-abi") == 0) {
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
			ret |= do_trans(infp, line, outfp, format);
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
		ret |= do_trans(infp, *argv, outfp, format);
		mfclose(infp);
	    }
	}
    } else {
	ret = do_trans(infp, "<stdin>", outfp, format);
    }

    return ret;
}

