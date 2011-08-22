/*
 * ======================================================================
 * This software has been created by Genome Research Limited (GRL).
 *
 * GRL hereby grants permission to use, copy, modify and distribute
 * this software and its documentation for non-commercial purposes
 * without fee at the user's own risk on the basis set out below.
 *
 * GRL neither undertakes nor accepts any duty whether contractual or
 * otherwise in connection with the software, its use or the use of
 * any derivative, and makes no representations or warranties, express
 * or implied, concerning the software, its suitability, fitness for
 * a particular purpose or non-infringement.
 *
 * In no event shall the authors of the software or GRL be responsible
 * or liable for any loss or damage whatsoever arising in any way
 * directly or indirectly out of the use of this software or its
 * derivatives, even if advised of the possibility of such damage.
 *
 * Our software can be freely distributed under the conditions set out
 * above, and must contain this copyright notice.
 * ======================================================================
 */

/*
 * This performs a linear (non-indexed) search for a trace in an SRF archive.
 *
 * It's not intended as a suitable production program or as a library of code
 * to use, but as a test and benchmark statistic.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <io_lib/Read.h>
#include <io_lib/misc.h>
#include <io_lib/ztr.h>
#include <io_lib/srf.h>

/* Fix MinGW's mkdir which only accepts one argument */
#if defined(__MINGW32__)
#    define mkdir(filename,mode) mkdir((filename))
#endif

#define CR 13            /* Decimal code of Carriage Return char */
#define LF 10            /* Decimal code of Line Feed char */
#define EOF_MARKER 26    /* Decimal code of DOS end-of-file marker */
#define MAX_REC_LEN 1024 /* Maximum size of input buffer */

#define SEQ (1 << 0)
#define PRB (1 << 1)
#define SIG2 (1 << 2)
#define INT (1 << 3)
#define NSE (1 << 4)
#define ALL (SEQ | PRB | SIG2 | INT | NSE)

#define SEQ_STR "seq"
#define PRB_STR "prb"
#define SIG2_STR "sig2"
#define INT_STR "int"
#define NSE_STR "nse"
#define ALL_STR "all"

#define TEXT (1 << 0)
#define SOLEXA (1 << 1)

#define TEXT_STR "text"
#define SOLEXA_STR "solexa"

#define CONSOLE_DEST (1 << 0)
#define FILE_DEST (1 << 1)
#define NONE_DEST (1 << 2)

#define CONSOLE_STR "console"
#define FILE_STR "file"
#define NONE_STR "none"

#define MAX_CYCLES 500

/* ------------------------------------------------------------------------ */

/*
 * All the reads and prefixes will be used to filter reads in an srf archive.
 * For reads the match with the reads in the archive must be exact.  Prefixes
 * need only match to the beginning of the read, so they are a way to filter
 * several reads that, for example, belong to the same lane and tile.  For
 * example, if the reads are <Center>:<Run>:<Lane>:<Tile>:<X>:<Y>, then a
 * prefix can progress from center all the way to y coordinate in terms of
 * read specificity.
 */
typedef struct {
    int prefixes_size;
    int reads_size;
    char** prefixes;  /* Prefixes to filter on. */
    char** reads;     /* Reads to filter on. */
} read_filter_t;

/* ------------------------------------------------------------------------ */

/* Find the first non-digit character in the name and make it variable. */
int parse_int_from_name(const char *name, int *variable, size_t *i, const char *variable_name) {
    size_t last = *i;
    int sep = 0;

    while (*i >= 0 && isdigit(name[*i])) {
	--(*i);
    }
    sep = name[*i];

    if(*i < 0 || last == *i) {
	fprintf(stderr, "Bad read name \"%s\" for \"%s\".  Read name needs to follow pattern of:\n", name, variable_name);
	fprintf(stderr, "\"prefix<separator><lane><separator><tile><separator><x><separator><y>\"\n");
	fprintf(stderr, "where lane, tile, x and y all consist of digits and the separators consist of non-numeric characters.\n");
	*variable = 0;
	return -1;
    }

    *variable = atoi(&name[(*i)+1]);
    --(*i);

    return sep;
}

/*
 * Parses a name assuming it consists of:
 * prefix<separator><lane><separator><tile><separator><x><separator><y>
 *
 * Sometimes though we have indexed samples, eg
 * "IL2_4381:1:1:1066:18864#43", where the last number has a different
 * meaning. We spot this by looking for a different separator.
 *
 * We fill out the supplied lane, tile, x and y parameters, or set them to
 * 0, 0, 0 and 0 if unknown.
 *
 * Returns 0 for success
 *        -1 for failure (unknown)
 */
int parse_name(char *name, int *lane, int *tile, int *x, int *y) {
    size_t len = strlen(name);
    size_t i = len-1;
    int sep1, sep2, sep3, sep4;

    /* Find the first non-digit character in the name and make it y. */
    if((sep1 = parse_int_from_name(name, y, &i, "y")) == -1)
	return -1;
    
    /* Find the second non-digit character in the name and make it x. */
    if((sep2 = parse_int_from_name(name, x, &i, "x")) == -1)
	return -1;

    /* Find the third non-digit character in the name and make it the tile. */
    if((sep3 = parse_int_from_name(name, tile, &i, "tile")) == -1)
	return -1;

    /* Find the fourth non-digit character in the name and make it the lane. */
    if((sep4 = parse_int_from_name(name, lane, &i, "lane")) == -1)
	return -1;

    /* Detect when the last value was something different, eg index */
    if (sep1 != sep2 && sep2 == sep3) {
	*y    = *x;
	*x    = *tile;
	*tile = *lane;
	sep1 = sep2;
	sep2 = sep3;
	sep3 = sep4;
	if ((sep4 = parse_int_from_name(name, lane, &i, "lane") == -1))
	    return -1;
    }

    if (sep1 != sep2 || sep2 != sep3) {
	static int done = 0;
	if (!done) {
	    done = 1;
	    fprintf(stderr, "Error: Name format unrecognised. "
		    "lane/tile/x/y values maybe incorrect.\n");
	}
    }

    return 0;
}

/*
 * Prints to the given file 4 values for every base in solexa format.  This is
 * used to output int, nse and sig2 type reads.
 */
void dump_samples4_solexa(FILE *fp, int lane, int tile, int x, int y,
			  int baseline, unsigned char *bytes, int nbytes) {
    int i, j, nb;
    int A[MAX_CYCLES], C[MAX_CYCLES], G[MAX_CYCLES], T[MAX_CYCLES];

    nb = nbytes/8;

    for (i = j = 0; i < nb; i++, j+= 2)
	A[i] = (bytes[j] << 8) | bytes[j+1];

    for (i = 0; i < nb; i++, j+= 2)
	C[i] = (bytes[j] << 8) | bytes[j+1];

    for (i = 0; i < nb; i++, j+= 2)
	G[i] = (bytes[j] << 8) | bytes[j+1];

    for (i = 0; i < nb; i++, j+= 2)
	T[i] = (bytes[j] << 8) | bytes[j+1];

    fprintf(fp, "%d\t%d\t%d\t%d", lane, tile, x, y);
    for (i = 0; i < nb; i++)
	fprintf(fp, "\t%d %d %d %d",
		A[i] - baseline,
		C[i] - baseline,
		G[i] - baseline,
		T[i] - baseline);
    fprintf(fp, "\n");
}

/*
 * Prints to the given file 4 values for every base in solexa format.  This is
 * used to output prb type reads.
 */
void dump_conf4_solexa(FILE *fp, char *seq, signed char *bytes, int nbytes) {
    int i, j, nb;
    int A[MAX_CYCLES], C[MAX_CYCLES], G[MAX_CYCLES], T[MAX_CYCLES];

    nb = nbytes/4;

    for (i = 0, j = nb; i < nb; i++) {
	switch(seq[i]) {
	case 'A': case 'a':
	    A[i] = bytes[i];
	    C[i] = bytes[j++];
	    G[i] = bytes[j++];
	    T[i] = bytes[j++];
	    break;
	case 'C': case 'c':
	    A[i] = bytes[j++];
	    C[i] = bytes[i];
	    G[i] = bytes[j++];
	    T[i] = bytes[j++];
	    break;
	case 'G': case 'g':
	    A[i] = bytes[j++];
	    C[i] = bytes[j++];
	    G[i] = bytes[i];
	    T[i] = bytes[j++];
	    break;
	default:
	    A[i] = bytes[j++];
	    C[i] = bytes[j++];
	    G[i] = bytes[j++];
	    T[i] = bytes[i];
	    break;
	}
    }

    for (i = 0; i < nb; i++) {
	if (i)
	    fputc('\t', fp);
	fprintf(fp, "%4d %4d %4d %4d", A[i], C[i], G[i], T[i]);
    }

    fprintf(fp, "\n");
}

/*
 * Prints a read in solexa format.  Depending on the given chunk type mode,
 * only some of the chunks are printed out for every read.  The files are
 * given in the following order: seq, prb, sig2, int, nse.  The files must
 * already be open.
 */
void dump_solexa(ztr_t *z, char *name, char mode, FILE **files) {
    int i, nc;
    ztr_chunk_t **chunks;
    char *seq;
    int lane, tile, x, y;
    parse_name(name, &lane, &tile, &x, &y);

    uncompress_ztr(z);

    chunks = ztr_find_chunks(z, ZTR_TYPE_BASE, &nc);
    if (nc != 1) {
	fprintf(stderr, "Zero or greater than one BASE chunks found.\n");
	return;
    }
    seq = chunks[0]->data+1;

    /* Sequence */
    if (mode & SEQ) {
	fprintf(files[0], "%d\t%d\t%d\t%d\t%.*s\n",
		lane, tile, x, y,
		chunks[0]->dlength-1,
		chunks[0]->data+1);
    }

    /* Confidence */
    if (mode & PRB) {
	chunks = ztr_find_chunks(z, ZTR_TYPE_CNF4, &nc);
	if (nc != 1) {
	    fprintf(stderr, "Zero or greater than one CNF chunks found.\n");
	    return;
	}

	dump_conf4_solexa(files[1], seq, chunks[0]->data+1, chunks[0]->dlength-1);
    }

    /* Traces */
    if (mode & SIG2) {
	chunks = ztr_find_chunks(z, ZTR_TYPE_SMP4, &nc);
	for (i = 0; i < nc; i++) {
	    char *key = ztr_lookup_mdata_value(z, chunks[i], "TYPE");
	    if (!key || 0 == strcmp(key, "PROC")) {
		key = ztr_lookup_mdata_value(z, chunks[i], "OFFS");
		dump_samples4_solexa(files[2], lane, tile, x, y, key ? atoi(key) : 0,
				     chunks[i]->data+2, chunks[i]->dlength-2);
		break;
	    }
	}
    }

    if (mode & INT) {
	chunks = ztr_find_chunks(z, ZTR_TYPE_SMP4, &nc);
	for (i = 0; i < nc; i++) {
	    char *key = ztr_lookup_mdata_value(z, chunks[i], "TYPE");
	    if (key && 0 == strcmp(key, "SLXI")) {
		key = ztr_lookup_mdata_value(z, chunks[i], "OFFS");
		dump_samples4_solexa(files[3], lane, tile, x, y, key ? atoi(key) : 0,
				     chunks[i]->data+2, chunks[i]->dlength-2);
		break;
	    }
	}
    }

    if (mode & NSE) {
	chunks = ztr_find_chunks(z, ZTR_TYPE_SMP4, &nc);
	for (i = 0; i < nc; i++) {
	    char *key = ztr_lookup_mdata_value(z, chunks[i], "TYPE");
	    if (key && 0 == strcmp(key, "SLXN")) {
		key = ztr_lookup_mdata_value(z, chunks[i], "OFFS");
		dump_samples4_solexa(files[4], lane, tile, x, y, key ? atoi(key) : 0,
				     chunks[i]->data+2, chunks[i]->dlength-2);
		break;
	    }
	}
    }
}

/*
 * Ripped out of io_lib's trace_dump program.
 * It reformats a trace to as printable ASCII.
 */
void dump_text(ztr_t *z, char *name, char mode, FILE **files) {
    Read *read;
    int i;

    uncompress_ztr(z);
    read = ztr2read(z); /* Inefficient; can do direct */

    if (read == NULL) {
        fprintf(stderr, "Tracedump was unable to open file %s\n", name );
        return;
    }

    fprintf(files[0], "[Trace]\n");
    fprintf(files[0], "%s\n", name);

    fprintf(files[0], "\n[Header]\n");
    fprintf(files[0], "%d\t\t# format\n",          read->format);
    fprintf(files[0], "%d\t\t# NPoints\n",         read->NPoints);
    fprintf(files[0], "%d\t\t# NBases\n",          read->NBases);
    fprintf(files[0], "%d\t\t# NFlows\n",          read->nflows);
    fprintf(files[0], "%d\t\t# maxTraceVal\n",     (int)read->maxTraceVal-read->baseline);
    fprintf(files[0], "%d\t\t# baseline\n",        read->baseline);
    fprintf(files[0], "%d\t\t# leftCutoff\n",      read->leftCutoff);
    fprintf(files[0], "%d\t\t# rightCutoff\n",     read->rightCutoff);

    fputs("\n[Bases]\n", files[0]);
    for (i = 0; i < read->NBases; i++) {
        fprintf(files[0], "%c %05d %+03d %+03d %+03d %+03d #%3d\n",
		read->base[i],
		read->basePos ? read->basePos[i] : 0,
		(int)read->prob_A[i],
		(int)read->prob_C[i],
		(int)read->prob_G[i],
		(int)read->prob_T[i],
		i);
    }

    if (read->NPoints) {
	fputs("\n[A_Trace]\n", files[0]);
	for(i = 0; i < read->NPoints; i++)
	    fprintf(files[0], "%d\t#%5d\n", (int)read->traceA[i] - read->baseline, i);

	fputs("\n[C_Trace]\n", files[0]);
	for(i = 0; i < read->NPoints; i++)
	    fprintf(files[0], "%d\t#%5d\n", (int)read->traceC[i] - read->baseline, i);

	fputs("\n[G_Trace]\n", files[0]);
	for(i = 0; i < read->NPoints; i++)
	    fprintf(files[0], "%d\t#%5d\n", (int)read->traceG[i] - read->baseline, i);

	fputs("\n[T_Trace]\n", files[0]);
	for(i = 0; i < read->NPoints; i++)
	    fprintf(files[0], "%d\t#%5d\n", (int)read->traceT[i] - read->baseline, i);
    }

    if (read->flow_order) {
        fputs("\n[Flows]\n", files[0]);
        for (i = 0; i < read->nflows; i++) {
            fprintf(files[0], "%c %5.2f  %u\t#%5d\n",
		    read->flow_order[i],
		    read->flow ? read->flow[i] : 0,
		    read->flow_raw ? read->flow_raw[i] : 0,
		    i);
        }
    }

    if (read->info) {
        fputs("\n[Info]\n", files[0]);
        fprintf(files[0], "%s\n", read->info);
    }

    read_deallocate(read);
}

/*
 * Print usage message to stderr and exit with the given \"code\".
 */
void usage(int code) {
    fprintf(stderr, "Usage: srf_dump_all [-c chunk_types] [-d destination_types]  [-f read_filter] [-n] [-o] [-t type_of_output] archive_name\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -c chunk_types\n");
    fprintf(stderr, "              Chunk types to output given as a comma delimited list of types.\n");
    fprintf(stderr, "              The following types are allowed: \"all\", \"seq\", \"prb\",\n");
    fprintf(stderr, "              \"int\", \"nse\", \"sig2\".\n");
    fprintf(stderr, "              The default is \"all\".\n");
    fprintf(stderr, "\n    -d destination_types\n");
    fprintf(stderr, "              Distinations to output to.\n");
    fprintf(stderr, "              The following types are allowed: \"console\", \"file\", \"none\".\n");
    fprintf(stderr, "              The default is \"console\".  Console and file can be together, \n"
	    "              but none is exclusive.\n");
    fprintf(stderr, "\n    -f read_filter\n");
    fprintf(stderr, "              The filter to apply to reads in the archive.  If reads match the\n");
    fprintf(stderr, "              filter they are dumped.\n");
    fprintf(stderr, "              The filter takes the form of <name>=<value>, where <name> can be\n");
    fprintf(stderr, "              \"read\", \"prefix\", \"file\".\n");
    fprintf(stderr, "              If the <name> is \"read\" the value is represents the name of a\n");
    fprintf(stderr, "              read to dump.  Only the matching read will be dumped.\n");
    fprintf(stderr, "              If the <name> is \"prefix\" the value is represents the prefix of\n");
    fprintf(stderr, "              the name of the reads to dump.  Only the matching reads will be\n              dumped.\n");
    fprintf(stderr, "              If the <name> is \"file\" the value is a file name where any\n");
    fprintf(stderr, "              number of \"read\" and \"prefix\" name value pairs can be included,\n              one per line.\n");
    fprintf(stderr, "              The default is no filter, which means all reads are dumped.\n");
    fprintf(stderr, "\n    -n        Output the total number of reads at the bottom.\n");
    fprintf(stderr, "\n    -o        Output the total number of reads only.  Nothing else is output.\n");
    fprintf(stderr, "\n    -t type_of_output\n");
    fprintf(stderr, "              Type of output.  Only one value allowed.\n");
    fprintf(stderr, "              Currently \"text\" and \"solexa\" is supported.\n");
    fprintf(stderr, "              The default is \"text\".\n");
    fprintf(stderr, "\n    -v        Print verbose messages.\n");
    fprintf(stderr, "\n");

    exit(code);
}

/*
 * Reads the lines from the fiven \"input\" file and creates the reads and prefixes
 * for the read filter.  The \"read_filter\" given is allocated.
 *
 * Returns 0 on success.  Errors cause the usage message and exit. 
 */
int read_filter_from_file(FILE *input, read_filter_t *read_filter)
{
    int   isNewline;              /* Boolean indicating we've read a CR or LF */
    long  lFileLen;               /* Length of file */
    long  lIndex;                 /* Index into cThisLine array */
    long  lLineCount;             /* Current line number */
    long  lLineLen;               /* Current line length */
    long  lStartPos;              /* Offset of start of current line */
    long  lTotalChars;            /* Total characters read */
    char  cThisLine[MAX_REC_LEN]; /* Contents of current line */
    char *cFile;                  /* Dynamically allocated buffer (entire file) */
    char *cThisPtr;               /* Pointer to current position in cFile */

    char *filter_type;
    char *prefix;
    char *read;

    fseek(input, 0L, SEEK_END);  /* Position to end of file */
    lFileLen = ftell(input);     /* Get file length */
    rewind(input);               /* Back to start of file */

    cFile = calloc(lFileLen + 1, sizeof(char));

    if(cFile == NULL )
	{
	    fprintf(stderr, "\nInsufficient memory to read file.\n");
	    return 0;
	}

    fread(cFile, lFileLen, 1, input); /* Read the entire file into cFile */

    lLineCount  = 0L;
    lTotalChars = 0L;

    cThisPtr    = cFile;              /* Point to beginning of array */

    while (*cThisPtr)                 /* Read until reaching null char */
	{
	    lIndex    = 0L;                 /* Reset counters and flags */
	    isNewline = 0;
	    lStartPos = lTotalChars;

	    while (*cThisPtr)               /* Read until reaching null char */
		{
		    if (!isNewline)               /* Haven't read a CR or LF yet */
			{
			    if (*cThisPtr == CR || *cThisPtr == LF) /* This char IS a CR or LF */
				isNewline = 1;                        /* Set flag */
			}

		    else if (*cThisPtr != CR && *cThisPtr != LF) /* Already found CR or LF */
			break;                                     /* Done with line */

		    /* Don't copy LS or CR */
		    if (*cThisPtr != CR && *cThisPtr != LF) {
			cThisLine[lIndex++] = *cThisPtr++; /* Add char to output and increment */
			++lTotalChars;
		    } else {
			cThisPtr++;
		    }

		} /* end while (*cThisPtr) */

	    cThisLine[lIndex] = '\0';     /* Terminate the string */
	    ++lLineCount;                 /* Increment the line counter */
	    lLineLen = strlen(cThisLine); /* Get length of line */

	    /* Find the one and only = in the string. */
	    if(strchr(cThisLine,'=') != NULL && (strchr(cThisLine,'=') == strrchr(cThisLine,'='))) {
		filter_type = strtok (cThisLine, "=");
	    } else {
		fprintf(stderr, "Baddly formatted read filter \"%s\".  Expected an \"=\" character in middle of filter.\n", cThisLine);
		usage(1);
	    }

	    if (!strcmp(filter_type, "prefix")) {
		prefix = strtok (NULL, "=");
		if(prefix == NULL) {
		    fprintf(stderr, "Bad prefix \"%s\" in read filter \"%s\".\n", prefix, cThisLine);
		    usage(1);
		} else {
		    ++(read_filter->prefixes_size);
		    read_filter->prefixes = (char**) realloc (read_filter->prefixes, read_filter->prefixes_size * sizeof(char *));
		    read_filter->prefixes[read_filter->prefixes_size - 1] =  (char*) calloc (strlen(prefix) + 1,sizeof(char));
		    strcpy(read_filter->prefixes[read_filter->prefixes_size - 1], prefix);
		}
	    } else if (!strcmp(filter_type, "read")) {
		read = strtok (NULL, "=");
		if(read == NULL) {
		    fprintf(stderr, "Bad read \"%s\" in read filter \"%s\".\n", read, cThisLine);
		    usage(1);
		} else {
		    ++(read_filter->reads_size);
		    read_filter->reads = (char**) realloc (read_filter->reads, read_filter->reads_size * sizeof(char *));
		    read_filter->reads[read_filter->reads_size - 1] =  (char*) calloc (strlen(read) + 1,sizeof(char));
		    strcpy(read_filter->reads[read_filter->reads_size - 1], read);
		}
	    } else {
		fprintf(stderr, "Unrecognized filter type \"%s\" given as part of read filter \"%s\".  The valid filter types are \"%s\".\n", filter_type, cThisLine, "prefix or read");
		usage(1);
	    }

	}

    free(cFile);
    return 0;
}

/*
 * Parse the given \"filter_value\" string for the filter value.
 *
 * Returns a read filter allocated on the heap which needs to be freed by the
 * calling function. Errors usually cause the usage message and exit.
 */
read_filter_t *get_read_filter(char *filter_value)
{
    char *filter_type = NULL;
    char *file_name = NULL;
    FILE *fp = NULL;
    char *prefix = NULL;
    char *read = NULL;

    /* Create read filter. */
    read_filter_t *read_filter = (read_filter_t *)calloc(1, sizeof(read_filter_t));
    if(read_filter == NULL) {
	return NULL;
    }
    read_filter->prefixes_size = 0;
    read_filter->reads_size = 0;

    /* Find the one and only = in the string. */
    if( strchr(filter_value,'=') != NULL &&
	(strchr(filter_value,'=') == strrchr(filter_value,'=')) ) {
	filter_type = strtok (filter_value,"=");
    } else {
        fprintf(stderr, "Baddly formatted read filter \"%s\".  Expected an \"=\" character in middle of filter.\n", filter_value);
        usage(1);
    }

    /* Check the string before the = is a valid filter type. */
    if(!strcmp("file", filter_type)) {

        /* Read the file. */
	file_name = strtok (NULL, "=");
	if(file_name == NULL) {
	    fprintf(stderr, "Bad file name \"%s\" in read filter \"%s\".\n", file_name, filter_value);
	    usage(1);
	}
	fp = fopen(file_name, "r");
	if(fp == NULL) {
            fprintf(stderr, "Bad file name \"%s\" in read filter \"%s\".\n", file_name, filter_value);
	    usage(1);
	}

	/* Read line by line. */
	if(read_filter_from_file(fp, read_filter)) {
	    fprintf(stderr, "Bad contents of file %s.\n", file_name);
	    usage(1);
	}

    } else if (!strcmp("prefix", filter_type)) {
	prefix = strtok (NULL, "=");
	if(prefix == NULL) {
	    fprintf(stderr, "Bad prefix \"%s\" in read filter \"%s\".\n", prefix, filter_value);
	    usage(1);
	} else {
	    ++(read_filter->prefixes_size);
	    read_filter->prefixes = (char**) malloc (read_filter->prefixes_size * sizeof(char *));
	    read_filter->prefixes[read_filter->prefixes_size - 1] =  (char*) calloc (strlen(prefix) + 1,sizeof(char));
	    strcpy(read_filter->prefixes[read_filter->prefixes_size - 1], prefix);
	}
    } else if (!strcmp("read", filter_type)) {
	read = strtok (NULL, "=");
	if(read == NULL) {
	    fprintf(stderr, "Bad read \"%s\" in read filter \"%s\".\n", read, filter_value);
	    usage(1);
	} else {
	    ++(read_filter->reads_size);
	    read_filter->reads = (char**) malloc (read_filter->reads_size * sizeof(char *));
	    read_filter->reads[read_filter->reads_size - 1] =  (char*) calloc (strlen(read) + 1, sizeof(char));
	    strcpy(read_filter->reads[read_filter->reads_size - 1], read);
	}
    } else {
	fprintf(stderr, "Unrecognized filter type \"%s\" given as part of read filter \"%s\".  The valid filter types are \"%s\".\n", filter_type, filter_value, "prefix, read, or file");
	usage(1);
    }

    return read_filter;
}

/*
 * Parse the comma delimited list of chunk types and put them in the single character \"mode\".
 *
 * Returns 0 on success.
 */
int get_chunk_types(char *arg, char *mode) {
    int num_allowed_types = 6;
    char *allowed_str_types[] = {SEQ_STR, INT_STR, NSE_STR, PRB_STR, SIG2_STR, ALL_STR};
    char allowed_types[] = {SEQ, INT, NSE, PRB, SIG2, ALL};
    char *type;
    int i = 0;

    type = strtok (arg,",");
    while(type) {
	for(i = 0; i < num_allowed_types; i++) {
	    if(!strcmp(type, allowed_str_types[i]) && !(*mode & allowed_types[i])) {
		*mode += allowed_types[i];
		break;
	    }
	}
        type = strtok (NULL, ",");
    }

    return 0;
}

/*
 * Parse the output format type and put it in the single character \"mode\".
 *
 * Returns 0 on success.
 */
int get_type_of_output(char *arg, char *mode) {
    int num_allowed_types = 2;
    char *allowed_str_types[] = {TEXT_STR, SOLEXA_STR};
    char allowed_types[] = {TEXT, SOLEXA};
    int i = 0;

    for(i = 0; i < num_allowed_types; i++) {
	if(!strcmp(arg, allowed_str_types[i])) {
	    *mode = allowed_types[i];
	    break;
	}
    }

    return 0;
}

/*
 * Parse the comma delimited list of destination types and put them in the single characted "\mode\".
 *
 * Returns 0 on success.
 */
int get_destination_types(char *arg, char *mode) {
    int num_allowed_types = 6;
    char *allowed_str_types[] = {CONSOLE_STR, FILE_STR, NONE_STR};
    char allowed_types[] = {CONSOLE_DEST, FILE_DEST, NONE_DEST};
    char *type;
    int i = 0;

    type = strtok (arg,",");
    while(type) {
	for(i = 0; i < num_allowed_types; i++) {
	    if(!strcmp(type, allowed_str_types[i]) && !(*mode & allowed_types[i])) {
		*mode += allowed_types[i];
		break;
	    }
	}
        type = strtok (NULL, ",");
    }

    /* NONE_DEST is exclusive, ignore others. */
    if(*mode & NONE_DEST) {
	*mode = NONE_DEST;
    }

    return 0;
}

/*
 * Returns 1 is the read \"name\" matches any of the reads or prefixes in the \"read_filter\".
 */
int check_read_name(read_filter_t *read_filter, char *name) {
    int i;

    for(i = 0; i < read_filter->prefixes_size; i++) {
	if(name == strstr(name, read_filter->prefixes[i])) {
	    return 1;
	}
    }

    for(i = 0; i < read_filter->reads_size; i++) {
	if(!strcmp(name, read_filter->reads[i])) {
	    free(read_filter->reads[i]);
	    read_filter->reads[i] = read_filter->reads[read_filter->reads_size - 1];
	    read_filter->reads_size--;
	    return 1;
	}
    }

    return 0;
}

void dump_read_filter(read_filter_t *read_filter) {
    int i;

    printf("Read filter used:\n");

    for(i = 0; i < read_filter->prefixes_size; i++) {
	printf("\tprefix[%d] = %s\n", i, read_filter->prefixes[i]);
    }

    for(i = 0; i < read_filter->reads_size; i++) {
	printf("\tread[%d] = %s\n", i, read_filter->reads[i]);
    }
}

void dump_chunk_mode(char mode) {
    printf("mode: %d.\n", mode);

    if(mode & SEQ) {
	printf("SEQ required.\n");
    }
    if(mode & INT) {
	printf("INT required.\n");
    }
    if(mode & NSE) {
	printf("NSE required.\n");
    }
    if(mode & PRB) {
	printf("PRB required.\n");
    }
    if(mode & SIG2) {
	printf("SIG2 required.\n");
    }
}

void dump_type_mode(char mode) {
    printf("mode: %d.\n", mode);

    if(mode & TEXT) {
	printf("TEXT required.\n");
    }
    if(mode & SOLEXA) {
	printf("SOLEXA required.\n");
    }
}

void dump_destination_mode(char mode) {
    printf("mode: %d.\n", mode);

    if(mode & CONSOLE_DEST) {
	printf("CONSOLE required.\n");
    }
    if(mode & FILE_DEST) {
	printf("FILE required.\n");
    }
    if(mode & NONE_DEST) {
	printf("NONE required.\n");
    }
}

/*
 * Open the solexa format file given the directory, lane, tile, and file type.
 * Uses solexa file naming conventions to create the file.
 *
 * Returns the successfully opened file or NULL.
 */
FILE *fopen_slx(char *dir, char *type, int lane, int tile) {
    char fn[1024];
    FILE *fp;

    sprintf(fn, "%s/s_%d_%04d_%s.txt", dir, lane, tile, type);
    printf("Opening file %s\n", fn);
    if (NULL == (fp = fopen(fn, "w+"))) {
	perror(fn);
	return NULL;
    }

    return fp;
}

/*
 * Given the archive name (ar_name), the chunk types to output (chunk_mode),
 * and some other parameters such as the read filter, open various files and
 * dumps solexa formated output to those files.  The number of reads is updated
 * in the \"num_reads\" parameter.
 *
 * Returns 0 on success.
 */
int process_srf_to_solexa_files(char *ar_name, char chunk_mode, int num_reads_only_mode, int filter_mode, read_filter_t *read_filter, long *num_reads) {
    srf_t *srf;
    char name[1024], dir2[1024];
    ztr_t *ztr;
    int last_lane = 0, last_tile = 0;
    /*               seq   prb   sig2  int   nse */
    FILE *files[] = {NULL, NULL, NULL, NULL, NULL};

    if (NULL == (srf = srf_open(ar_name, "rb"))) {
	perror(ar_name);
	return 1;
    }

    char *cp = strrchr(ar_name, '.');
    if (cp) *cp = 0;
    sprintf(dir2, "%s.run", ar_name);
    mkdir(dir2, 0777);

    while (NULL != (ztr = srf_next_ztr(srf, name, 0))) {
	int lane, tile, x, y;

	parse_name(name, &lane, &tile, &x, &y);

	if (last_lane != lane || last_tile != tile) {
	    fprintf(stderr, "New tile: %d/%d\n", lane, tile);
	    last_lane = lane;
	    last_tile = tile;

	    if (files[0]) { fclose(files[0]); files[0] = NULL; }
	    if (files[1]) { fclose(files[1]); files[1] = NULL; }
	    if (files[2]) { fclose(files[2]); files[2] = NULL; }
	    if (files[3]) { fclose(files[3]); files[3] = NULL; }
	    if (files[4]) { fclose(files[4]); files[4] = NULL; }

	    if (chunk_mode & SEQ) {
		files[0] = fopen_slx(dir2, "seq", lane, tile);
	    }
	    if (chunk_mode & PRB) {
		files[1] = fopen_slx(dir2, "prb", lane, tile);
	    }
	    if (chunk_mode & SIG2) {
		files[2] = fopen_slx(dir2, "sig2", lane, tile);
	    }
	    if (chunk_mode & INT) {
		files[3] = fopen_slx(dir2, "int", lane, tile);
	    }
	    if (chunk_mode & NSE) {
		files[4] = fopen_slx(dir2, "nse", lane, tile);
	    }
	}

	if(!num_reads_only_mode) {
	    if(!filter_mode || (filter_mode && check_read_name(read_filter, name))) {
		dump_solexa(ztr, name, chunk_mode, files);
		++(*num_reads);
	    }
	} else {
	    ++(*num_reads);
	}

	delete_ztr(ztr);
    }

    if (files[0]) fclose(files[0]);
    if (files[1]) fclose(files[1]);
    if (files[2]) fclose(files[2]);
    if (files[3]) fclose(files[3]);
    if (files[4]) fclose(files[4]);

    srf_destroy(srf, 1);
    return 0;
}

/*
 * Open the text format file given the directory, lane and tile.
 * Uses solexa file naming conventions to create the file, except that the file
 * has a \"dump.txt\" type.
 *
 * Returns the successfully opened file or NULL.
 */
FILE *fopen_text(char *dir, int lane, int tile) {
    char fn[1024];
    FILE *fp;

    sprintf(fn, "%s/s_%d_%04d_dump.txt", dir, lane, tile);
    if (NULL == (fp = fopen(fn, "w+"))) {
	perror(fn);
	return NULL;
    }

    return fp;
}

/*
 * Given the archive name (ar_name), the chunk types to output (chunk_mode),
 * and some other parameters such as the read filter, open one file and dumps
 * text formatted output to that file.  The number of reads is updated in the
 * \"num_reads\" parameter.
 *
 * Returns 0 on success.
 */
int process_srf_to_text_files(char *ar_name, char chunk_mode, int num_reads_only_mode, int filter_mode, read_filter_t *read_filter, long *num_reads) {
    srf_t *srf;
    char name[1024], dir2[1024];
    ztr_t *ztr;
    int last_lane = 0, last_tile = 0;
    /*               dump */
    FILE *files[] = {NULL};

    if (NULL == (srf = srf_open(ar_name, "rb"))) {
	perror(ar_name);
	return 1;
    }

    char *cp = strrchr(ar_name, '.');
    if (cp) *cp = 0;
    sprintf(dir2, "%s.run", ar_name);
    mkdir(dir2, 0777);

    while (NULL != (ztr = srf_next_ztr(srf, name, 0))) {
	int lane, tile, x, y;

	parse_name(name, &lane, &tile, &x, &y);

	if (last_lane != lane || last_tile != tile) {
	    fprintf(stderr, "New tile: %d/%d\n", lane, tile);
	    last_lane = lane;
	    last_tile = tile;

	    if (files[0]) { fclose(files[0]); files[0] = NULL; }

	    if (chunk_mode) {
		files[0] = fopen_text(dir2, lane, tile);
	    }
	}

	if(!num_reads_only_mode) {
	    if(!filter_mode || (filter_mode && check_read_name(read_filter, name))) {
		dump_text(ztr, name, chunk_mode, files);
		++(*num_reads);
	    }
	} else {
	    ++(*num_reads);
	}

	delete_ztr(ztr);
    }

    if (files[0]) fclose(files[0]);

    srf_destroy(srf, 1);
    return 0;
}

/* ------------------------------------------------------------------------ */

/*
 * Main method.
 */
int main(int argc, char **argv) {
    int i;
    long num_reads = 0;
    int num_reads_mode = 0;
    int num_reads_only_mode = 0;
    int filter_mode = 0;
    char *ar_name = NULL;
    srf_t *srf = NULL;
    char name[512];
    ztr_t *ztr = NULL;
    read_filter_t *read_filter = NULL;
    char *filter_value = NULL;

    int c;
    int errflg = 0;
    extern char *optarg;
    extern int optind, optopt;

    char chunk_mode = ALL;
    char type_mode = TEXT;
    char destination_mode = CONSOLE_DEST;
    FILE **files = NULL;
    int verbose = 0;

    if (argc < 2) {
	fprintf(stderr, "Please specify an archive name.\n");
	usage(1);
    }

    while ((c = getopt(argc, argv, ":c:d:f:not:v")) != -1) {
        switch (c) {
        case 'c':
	    chunk_mode = 0;
	    if(get_chunk_types(optarg, &chunk_mode) || !chunk_mode) {
                fprintf(stderr,
			"Invalid value \"%s\" given to option -%c.\n", optarg, c);
		errflg++;
	    }
	    break;
        case 'd':
	    destination_mode = 0;
	    if(get_destination_types(optarg, &destination_mode) || !destination_mode) {
                fprintf(stderr,
			"Invalid value \"%s\" given to option -%c.\n", optarg, c);
		errflg++;
	    }
	    break;
        case 'f':
	    if (num_reads_only_mode) {
		fprintf(stderr,
			"Option -%c is exclusing with option -o.\n", c);
                errflg++;
	    } else {
		filter_mode++;
		filter_value = optarg;
	    }
            break;
        case 'n':
	    if (num_reads_only_mode) {
		fprintf(stderr,
			"Option -%c is exclusing with option -o.\n", c);
                errflg++;
	    }
            else
                num_reads_mode++;
            break;
        case 'o':
	    if (num_reads_mode) {
		fprintf(stderr,
			"Option -%c is exclusing with option -n.\n", c);
                errflg++;
	    } else if (filter_mode) {
		fprintf(stderr,
			"Option -%c is exclusing with option -f.\n", c);
                errflg++;
	    } else
                num_reads_only_mode++;
            break;
        case 't':
	    type_mode=0;
	    if(get_type_of_output(optarg, &type_mode) || !type_mode) {
                fprintf(stderr,
			"Invalid value \"%s\" given to option -%c.\n", optarg, c);
		errflg++;
	    }
	    break;
        case 'v':
	    verbose++;
            break;
        case ':':       /* -? without operand */
            fprintf(stderr,
                    "Option -%c requires an operand\n", optopt);
            errflg++;
            break;
        case '?':
            fprintf(stderr,
                    "Unrecognised option: -%c\n", optopt);
            errflg++;
        }
    }

    if (errflg) {
	usage(1);
    }

    if(optind + 1 != argc) {
        fprintf(stderr, "The archive name must be the last argument.\n");
        usage(1);
    } else {
        ar_name = argv[optind];
    }
    
    if(filter_mode) {
	read_filter = get_read_filter(filter_value);
	if(verbose) {
	    dump_read_filter(read_filter);
	}
    }

    if(chunk_mode && verbose) {
	dump_chunk_mode(chunk_mode);
    }

    if(type_mode && verbose) {
	dump_type_mode(type_mode);
    }

    if(destination_mode && verbose) {
	dump_destination_mode(destination_mode);
    }

    if (!access(ar_name, R_OK)) {
	if(verbose) {
	    printf("Dumping from archive %s.\n", ar_name);
	}
    } else {
        fprintf(stderr, "Archive %s not found.\n", ar_name);
    }

    if(destination_mode & NONE_DEST) {
	if (NULL == (srf = srf_open(ar_name, "rb"))) {
	    perror(ar_name);
	    return 1;
	}

	while (NULL != (ztr = srf_next_ztr(srf, name, 0))) {
	    if(!num_reads_only_mode) {
		if(!filter_mode || (filter_mode && check_read_name(read_filter, name))) {
		    ++num_reads;
		}
	    } else {
		++num_reads;
	    }
	    delete_ztr(ztr);
	}
    }

    if(destination_mode & CONSOLE_DEST) {
	/* Cosole means stdout. */
	files = malloc((sizeof(FILE *))*5);
	for(i = 0; i < 5; i++) {
	    files[i] = stdout;
	}

	if (NULL == (srf = srf_open(ar_name, "rb"))) {
	    perror(ar_name);
	    return 1;
	}

	while (NULL != (ztr = srf_next_ztr(srf, name, 0))) {
	    if(!num_reads_only_mode) {
		if(!filter_mode || (filter_mode && check_read_name(read_filter, name))) {
		    if(type_mode & SOLEXA) {
			dump_solexa(ztr, name, chunk_mode, files);
		    } else if(type_mode & TEXT) {
			dump_text(ztr, name, chunk_mode, files);
		    } else {
			fprintf(stderr, "Assertion error on type_mode (%c).\nExiting.\n", type_mode);
			exit(1);
		    }
		    ++num_reads;
		}
	    } else {
		++num_reads;
	    }
	    delete_ztr(ztr);
	}
    }

    if(destination_mode & FILE_DEST) {
	if(type_mode & SOLEXA) {
	    process_srf_to_solexa_files(ar_name, chunk_mode, num_reads_only_mode, filter_mode, read_filter, &num_reads);
	} else if(type_mode & TEXT) {
	    process_srf_to_text_files(ar_name, chunk_mode, num_reads_only_mode, filter_mode, read_filter, &num_reads);
	} else {
	    fprintf(stderr, "Assertion error on type_mode (%c).\nExiting.\n", type_mode);
	    exit(1);
	}
    }

    if(num_reads_mode || num_reads_only_mode) {
	printf("\nReads: %ld\n", num_reads);
    }

    srf_destroy(srf, 1);

    return 0;
}
