#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>

#include <io_lib/Read.h>
#include <io_lib/traceType.h>
#include <io_lib/seqIOABI.h>
#include <io_lib/open_trace_file.h>

static char const rcsid[] = "$Id: convert_trace.c,v 1.12 2008-02-20 16:07:44 jkbonfield Exp $";

#define MAX(a,b) ((a)>(b)?(a):(b))

struct opts {
    char *name;
    char *fofn;
    char *passed;
    char *failed;
    char *error;
    int in_format;
    int out_format;
    int scale;
    int sub_background;
    int subtract;
    int normalise;
    int min_normalise;
    int compress_mode;
    int dots;
    int noneg;
    int signed_trace;
    int skipx;
    int start;
    int end;
};

/*
 * Removes any negative values from a trace by moving each channel
 * independently so that its lowest value is 0.
 */
void noneg(Read *r) {
    int i, j, k;
    signed int min;
    TRACE *t;

    /* Find the real end of the data */
    for (k = r->NPoints-1; k >= 0; k--)
	if (r->traceA[k] ||
	    r->traceC[k] ||
	    r->traceG[k] ||
	    r->traceT[k])
	    break;

    for (j = 0; j < 4; j++) {
	switch(j) {
	case 0:
	    t = r->traceA;
	    break;
	case 1:
	    t = r->traceC;
	    break;
	case 2:
	    t = r->traceG;
	    break;
	case 3:
	default:
	    t = r->traceT;
	    break;
	}

	/* Find the lowest -ve value per lane */
	for (min = i = 0; i <= k; i++) {
	    if (min > ((int16_t *)t)[i])
		min = ((int16_t *)t)[i];
	}

	/* And shift everything back up */
	for (i = 0; i <= k; i++) {
	    t[i] -= min;
	}
    }
}


/*
 * Removes any negative values from a trace by moving the trace up so that
 * the lowest overall value is 0. This differs to noneg above by using
 * a global shift for all channels and also setting the read 'baseline'.
 */
void signed_trace(Read *r) {
    int i, k;
    signed int min;

    /* Find the real end of the data */
    for (k = r->NPoints-1; k >= 0; k--)
	if (r->traceA[k] ||
	    r->traceC[k] ||
	    r->traceG[k] ||
	    r->traceT[k])
	    break;

    /* Find the lowest -ve value per lane */
    for (min = i = 0; i <= k; i++) {
	if (min > ((int16_t *)(r->traceA))[i])
	    min = ((int16_t *)(r->traceA))[i];
	if (min > ((int16_t *)(r->traceC))[i])
	    min = ((int16_t *)(r->traceC))[i];
	if (min > ((int16_t *)(r->traceG))[i])
	    min = ((int16_t *)(r->traceG))[i];
	if (min > ((int16_t *)(r->traceT))[i])
	    min = ((int16_t *)(r->traceT))[i];
    }

    r->baseline = -min;

    /* And shift everything back up */
    for (i = 0; i <= k; i++) {
	r->traceA[i] -= min;
	r->traceC[i] -= min;
	r->traceG[i] -= min;
	r->traceT[i] -= min;
    }
}

/*
 * Scales trace values from 0 to scale, but only if they are larger.
 */
void rescale_trace(Read *r, int scale) {
    double s;
    int i;

    if (r->maxTraceVal <= scale)
	return;
    
    s = ((double)scale)/r->maxTraceVal;

    for (i = 0; i < r->NPoints; i++) {
	r->traceA[i] = r->traceA[i] * s + 0.5;
	r->traceC[i] = r->traceC[i] * s + 0.5;
	r->traceG[i] = r->traceG[i] * s + 0.5;
	r->traceT[i] = r->traceT[i] * s + 0.5;
    }

    r->maxTraceVal = scale;
}

#if 0
/* OLD method, treats all channels together and assumes the same baseline for
 * each
 */
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
#endif

static void subtract_background_ch(TRACE *channel, int nchannel) {
    int i, j, bg;
    int win_len = 501, win_len2 = win_len/2;
    TRACE *copy;

    if (NULL == (copy = (TRACE *)malloc(sizeof(*copy) * nchannel)))
	return;

    if (nchannel < win_len)
	win_len = nchannel;

    /* Take lowest background over win_len and subtract it */
    for (i = 0; i < nchannel; i++) {
	/* Could optimise this considerably */
	bg = INT_MAX;
	for (j = -win_len2; j < win_len2; j++) {
	    if (i+j < 0) continue;
	    if (i+j >= nchannel) break;
		
	    if (channel[i + j] < bg)
		bg = channel[i + j];
	}
	
	copy[i] = channel[i] - bg;
    }

    memcpy(channel, copy, nchannel * sizeof(*copy));
    free(copy);
}

/*
 * Find the average background level of a trace, and subtract this from the
 * peak heights.
 */
void subtract_background(Read *r) {
    subtract_background_ch(r->traceA, r->NPoints);
    subtract_background_ch(r->traceC, r->NPoints);
    subtract_background_ch(r->traceG, r->NPoints);
    subtract_background_ch(r->traceT, r->NPoints);
}

int int_compar(const void *a, const void *b) {
    return *(const TRACE *)a - *(const TRACE *)b;
}

int find_bg(TRACE *data, int ndata) {
    int i, bg;
    TRACE *copy = (TRACE *)malloc(ndata * sizeof(TRACE));

    /* Sort the trace samples by amplitude */
    memcpy(copy, data, ndata * sizeof(TRACE));
    qsort(copy, ndata, sizeof(TRACE), int_compar);

    /* Find the first non-zero value */
    for (i = 0; i < ndata && !copy[i]; i++)
	;

    /*
     * Now take a slie 0.05 through the remainder of the array and set this
     * as our background.
     */
    bg = copy[(int)((ndata - i) * 0.05 + i)];

    free(copy);
    return bg;
}

void trace_freq(TRACE *data, int ndata) {
    int i, bg;
    bg = find_bg(data, ndata);

    for (i = 0; i < ndata; i++) {
	data[i] = MAX(data[i] - bg, 0);
    }
}

/*
 * Separates out the dyes using a deconvolution matrix.
 * The order of elements in the matrix is C A G T.
 * A test matrix for the 373. Taken from the BASS distribution.
 */
double matrix[5][4] = {
  { 0.002439782,        -0.0015053751,       0.00011857301,    2.8906948e-06},
  {-0.00075353298,       0.0032971052,      -0.006198165,      0.00014828549},
  { 0.00020249287,      -0.0017620348,       0.010530438,     -0.0020235507 },
  {-0.001144423,        -4.857673e-06,      -0.0018845701,     0.00395431   },
  {-0.12451385,          0.368916,          -2.928292,        -3.3142638    }
};
void separate_dyes(Read *r, double M[][4]) {
    int i, j;

    for (i = 0; i < r->NPoints; i++) {
	int C, A, G, T;
	double sep[4];

	C = r->traceC[i];
	A = r->traceA[i];
	G = r->traceG[i];
	T = r->traceT[i];

	for (j = 0; j < 4; j++)
	  sep[j] = C*M[0][j] + A*M[1][j] + G*M[2][j] + T*M[3][j] + M[4][j];

	for (j = 0; j < 4; j++)
	    sep[j] += 10;

	/* hack!
	   sep[0] += 0.1;
	   sep[1] += -0.4;
	   sep[2] += 2.9;
	   sep[3] += 3.2;
	*/

	r->traceC[i] = sep[0] < 0 ? 0 : 1000 * sep[0];
	r->traceA[i] = sep[1] < 0 ? 0 : 1000 * sep[1];
	r->traceG[i] = sep[2] < 0 ? 0 : 1000 * sep[2];
	r->traceT[i] = sep[3] < 0 ? 0 : 1000 * sep[3];
    }
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

/*
 * Rescales peak heights based on a moving "marker". The marker tracks
 * up and down (attack and decay) based on the difference between itself and
 * the trace envelope. We then divide by the marker value to attempt to
 * normalise peak heights.
 *
 * min_marker is used to avoid scaling up noise and represents the minimum
 * value the marker is allowed to reach. Make sure it is > 0 or divide by
 * zero may occur.
 */
void rescale_heights(Read *r, int min_marker) {
    double marker = 0;
    int i, j, max, mtv = 0;
    TRACE *tx[4];

    tx[0] = r->traceA;
    tx[1] = r->traceC;
    tx[2] = r->traceG;
    tx[3] = r->traceT;

    for (i = 0; i < r->NPoints; i++) {
	for (max = j = 0; j < 4; j++)
	    if (max < tx[j][i])
		max = tx[j][i];
	if (!marker) {
	    marker = max;
	} else {
	    if (max >= marker) {
		/* attack */
		marker += (max - marker) / 20.0;
	    } else {
		/* decay */
		marker -= (marker - max) / 10.0;
	    }
	}
	if (marker < min_marker)
	    marker = min_marker;

	for (j = 0; j < 4; j++) {
	    double new = tx[j][i] * 2000.0/marker;
	    tx[j][i] = new > 32767 ? 32767 : new;
	    if (mtv < tx[j][i])
		mtv = tx[j][i];
	}

    }

    r->maxTraceVal = mtv;
}

/* Removes every other sample as a crude way to reduce file size */
void skipx(Read *r) {
    int i, j;
    for (i = j = 0; j < r->NPoints/2; i+=2, j++) {
	r->traceA[j] = (r->traceA[i] + r->traceA[i+1]) / 2;
	r->traceC[j] = (r->traceC[i] + r->traceC[i+1]) / 2;
	r->traceG[j] = (r->traceG[i] + r->traceG[i+1]) / 2;
	r->traceT[j] = (r->traceT[i] + r->traceT[i+1]) / 2;
    }
    r->NPoints = j;

    for (i = 0; i < r->NBases; i++) {
	r->basePos[i] /= 2;
    }
}

void clip_range(Read *r, int left, int right) {
    int i, j;
    if (left != -1) {
	for (i = 0, j = left; j < r->NPoints; i++, j++) {
	    r->traceA[i] = r->traceA[j];
	    r->traceC[i] = r->traceC[j];
	    r->traceG[i] = r->traceG[j];
	    r->traceT[i] = r->traceT[j];
	}
	right -= left;
    }
    if (right > 0) {
	r->NPoints = right;
    }
}


int convert(mFILE *infp, mFILE *outfp, char *infname, char *outfname,
	    struct opts *opts) {
    Read *r;

    if (NULL == (r = mfread_reading(infp, infname, opts->in_format))) {
	fprintf(stderr, "failed to read file %s\n", infname);
	return 1;
    }

    if (opts->start != -1 || opts->end != -1)
	clip_range(r, opts->start, opts->end);

    if (opts->skipx) {
	skipx(r);
    }

    if (opts->noneg)
	noneg(r);

    if (opts->signed_trace)
	signed_trace(r);

    if (opts->subtract) {
	int i;
	for (i = 0; i < r->NPoints; i++) {
	    r->traceA[i] = MAX(0, r->traceA[i] - opts->subtract);
	    r->traceC[i] = MAX(0, r->traceC[i] - opts->subtract);
	    r->traceG[i] = MAX(0, r->traceG[i] - opts->subtract);
	    r->traceT[i] = MAX(0, r->traceT[i] - opts->subtract);
	}
    }

    if (opts->sub_background) { 
	/*
	trace_freq(r->traceA, r->NPoints);	
	trace_freq(r->traceC, r->NPoints);	
	trace_freq(r->traceG, r->NPoints);	
	trace_freq(r->traceT, r->NPoints);	
	*/
	subtract_background(r);
	/*
	separate_dyes(r, matrix);
	trace_freq(r->traceA, r->NPoints);	
	trace_freq(r->traceC, r->NPoints);	
	trace_freq(r->traceG, r->NPoints);	
	trace_freq(r->traceT, r->NPoints);
	*/
	reset_max_called_height(r);
    }

    if (opts->normalise) {
	rescale_heights(r, opts->min_normalise);
    }

    if (opts->scale) {
	rescale_trace(r, opts->scale);
    }

    if (opts->name)
	r->ident = strdup(opts->name);
    else if (0 == strcmp(outfname, "(stdout)"))
	r->ident = strdup(infname);
    else
	r->ident = strdup(outfname);

    if (opts->compress_mode != -1)
	set_compression_method(opts->compress_mode);

    if (0 != (mfwrite_reading(outfp, r, opts->out_format))) {
	fprintf(stderr, "failed to write file %s\n", outfname);
	read_deallocate(r);
	return 1;
    }

    read_deallocate(r);
    return 0;
}


void usage(void) {
    puts("Usage: convert_trace [options] [informat outformat] < in > out");
    puts("Or     convert_trace [options] -fofn file_of_filenames");
    puts("\nOptions are:");
    puts("    -in_format format         Format for input (defaults to any");
    puts("    -out_format format        Format for output (default ztr)");
    puts("    -fofn file_of_filenames   Get \"Input Output\" names from a fofn");
    puts("    -passed fofn              Output fofn of passed names");  
    puts("    -error errs               Redirect stderr to file \"errs\"");
    puts("    -failed fofn              Output fofn of failed names");  
    puts("    -name id                  ID line for experiment file output");
    puts("    -subtract_background      Auto-subtracts the trace background");
    puts("    -subtract amount          Subtracts a specified background amount");
    puts("    -normalise                Normalises peak heights");
    puts("    -min_normalise            Minimum trace amp for normalising");
    puts("    -scale range              Downscales peaks to 0-range");
    puts("    -compress mode            Compress file output (not if stdout)");
    puts("    -abi_data counts          ABI DATA lanes to copy: eg 9,10,11,12");
    puts("    -signed                   Apply global shift to avoid negative values");
    puts("    -noneg                    Shift each channel independently to avoid -ve");
    puts("    --                        Explicitly state end of options");
    exit(1);
}

int main(int argc, char **argv) {
    struct opts opts;

    opts.in_format = TT_ANY;
    opts.out_format = TT_ZTR;
    opts.scale = 0;
    opts.sub_background = 0;
    opts.subtract = 0;
    opts.normalise = 0;
    opts.min_normalise = 100;
    opts.name = NULL;
    opts.compress_mode = -1;
    opts.dots = 0;
    opts.noneg = 0;
    opts.signed_trace = 0;
    opts.fofn = NULL;
    opts.passed = NULL;
    opts.failed = NULL;
    opts.error = NULL;
    opts.skipx = 0;
    opts.start = -1;
    opts.end = -1;
    
    for (argc--, argv++; argc > 0; argc--, argv++) {
	if (**argv != '-')
	    break;

	if (strcmp(*argv, "-start") == 0) {
	    opts.start = atoi(*++argv);
	    argc--;

	} else if (strcmp(*argv, "-end") == 0) {
	    opts.end = atoi(*++argv);
	    argc--;

	} else if (strcmp(*argv, "-scale") == 0) {
	    opts.scale = atoi(*++argv);
	    argc--;

	} else if (strcmp(*argv, "-fofn") == 0) {
	    opts.fofn = *++argv;
	    argc--;

	} else if (strcmp(*argv, "-passed") == 0) {
	    opts.passed = *++argv;
	    argc--;

	} else if (strcmp(*argv, "-failed") == 0) {
	    opts.failed = *++argv;
	    argc--;

	} else if (strcmp(*argv, "-error") == 0) {
	    opts.error = *++argv;
	    argc--;

	} else if (strcmp(*argv, "-subtract_background") == 0) {
	    opts.sub_background = 1;

	} else if (strcmp(*argv, "-subtract") == 0) {
	    opts.subtract = atoi(*++argv);
	    argc--;

	} else if (strcmp(*argv, "-normalise") == 0) {
	    opts.normalise = 1;

	} else if (strcmp(*argv, "-min_normalise") == 0) {
	    opts.min_normalise = atoi(*++argv);
	    argc--;

	} else if (strcmp(*argv, "-dots") == 0) {
	    opts.dots = 1;

	} else if (strcmp(*argv, "-noneg") == 0) {
	    opts.noneg = 1;

	} else if (strcmp(*argv, "-signed") == 0) {
	    opts.signed_trace = 1;

	} else if (strcmp(*argv, "-skipx") == 0) {
	    opts.skipx = 1;

	} else if (strcmp(*argv, "-in_format") == 0) {
	    argv++;
	    argc--;
	    if (TT_UNK == (opts.in_format = trace_type_str2int(*argv)))
		opts.in_format = atoi(*argv);

	} else if (strcmp(*argv, "-name") == 0) {
	    opts.name = *++argv;
	    argc--;

	} else if (strcmp(*argv, "-out_format") == 0) {
	    argv++;
	    argc--;
	    if (TT_UNK == (opts.out_format = trace_type_str2int(*argv)))
		opts.out_format = atoi(*argv);

	} else if (strcmp(*argv, "-compress") == 0) {
	    opts.compress_mode = compress_str2int(*++argv);
	    argc--;

	} else if (strcmp(*argv, "-abi_data") == 0) {
	    int c1, c2, c3, c4;
	    argc--;
	    if (4 == sscanf(*++argv, "%d,%d,%d,%d", &c1, &c2, &c3, &c4)) {
		abi_set_data_counts(c1, c2, c3, c4);
	    } else {
		usage();
	    }

	} else if (strcmp(*argv, "--") == 0) {
	    break;

	} else {
	    usage();
	}
    }

    if (argc == 2) {
	/* Old syntax, for backwards compatibility */

	if (TT_UNK == (opts.in_format = trace_type_str2int(argv[0])))
	    opts.in_format = atoi(argv[0]);
	if (TT_UNK == (opts.out_format = trace_type_str2int(argv[1])))
	    opts.out_format = atoi(argv[1]);
    } else if (argc != 0) {
	usage();
    }


    /*
     * Added by SAK: Allow redirection of error output to file, due to
     * problems with Java exec
     */
    if (NULL != opts.error) {
	int fd;

	fprintf(stderr,"* Redirecting stderr to %s\n", opts.error);

	close(2); /* close fd with stderr */
	if (-1 == (fd = creat(opts.error, 0666))) {
	    exit(1);
	}
    }

    if (!opts.fofn) {
	return convert(mstdin(), mstdout(), "(stdin)", "(stdout)", &opts);
    }

    /* else */ {
	mFILE *fpin, *fpout;
	FILE *fppassed = NULL, *fpfailed = NULL;
	char *infname, *outfname;
	int ret, ret_all = 0;
	char line[8192], line2[8192];

	FILE *fofn_fp;

	if (NULL == (fofn_fp = fopen(opts.fofn, "r"))) {
	    perror(opts.fofn);
	    return -1;
	}

	if (opts.passed && NULL == (fppassed = fopen(opts.passed, "w"))) {
	    perror(opts.passed);
	    return -1;
	}

	if (opts.failed && NULL == (fpfailed = fopen(opts.failed, "w"))) {
	    perror(opts.failed);
	    return -1;
	}

	while (fgets(line, 8192, fofn_fp) != NULL) {
	    int i, j, len;
	    
	    /* Find input and output name, escaping spaces as needed */
	    len = strlen(line);
	    outfname = NULL;
	    for (i = j = 0; i < len; i++) {
		if (line[i] == '\\' && i != len-1) {
		    line2[j++] = line[++i];
		} else if (line[i] == ' ') {
		    line2[j++] = 0;
		    outfname = &line2[j];
		} else if (line[i] != '\n') {
		    line2[j++] = line[i];
		}
	    }
	    line2[j] = 0;
	    infname = line2;

	    /* Don't clobber input */
	    if (!strcmp(infname, outfname)) {
		fprintf(stderr,"* Inputfn %s == Outputfn %s ...skipping\n",
			infname, outfname);
		if (fpfailed)
		    fprintf(fpfailed, "%s\n", infname);
		continue;
	    }

	    /* Open input and output files */
	    if (opts.in_format == TT_EXP) {
		fpin = open_exp_mfile(infname, NULL);
	    } else {
		fpin = open_trace_mfile(infname, NULL);
	    }
	    if (NULL == fpin) {
		char buf[2048];
		sprintf(buf, "ERROR %s", infname);
		perror(buf);
		if (opts.dots) {
		    fputc('!', stdout);
		    fflush(stdout);
		}
		if (fpfailed)
		    fprintf(fpfailed, "%s\n", infname);
		continue;
	    }

	    if (outfname) {
		if (NULL == (fpout = mfopen(outfname, "wb+"))) {
		    char buf[2048];
		    sprintf(buf, "ERROR %s", outfname);
		    perror(buf);
		    mfclose(fpin);
		    if (opts.dots) {
			fputc('!', stdout);
			fflush(stdout);
		    }
		    if (fpfailed)
			fprintf(fpfailed, "%s\n", infname);
		    continue;
		}
	    } else {
		outfname = "(stdout)";
		fpout = mstdout();
	    }

	    /* Convert */
	    ret = convert(fpin, fpout, infname, outfname, &opts);
	    ret_all |= ret;
	    if (opts.dots) {
		fputc(ret ? '!' : '.', stdout);
		fflush(stdout);
	    }
	    if (ret) {
		if (fpfailed)
		    fprintf(fpfailed, "%s\n", infname);
	    } else {
		if (fppassed)
		    fprintf(fppassed, "%s\n", infname);
	    }

	    /* Tidy up */
	    mfclose(fpin);
	    if (fpout != mstdout())
		mfclose(fpout);
	}

	fclose(fofn_fp);

	return ret_all;
    }
}
