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
  Title:       seqIOPlain
  
  File: 	 seqIOPlain.c
  Purpose:	 IO of plain sequences
  Last update: Tuesday Jan 15 1991
  
  Change log:
  
  28.11.90 SD  put undesirables under STLOUIS compilation flag
  15.01.91 SD  new include file (opp.h)
  17.09.91 LFW changed STLOUIS compilation flag to SAVE_EDITS
  and AUTO_CLIP
  08.09.94 JKB Plain files now also uses the ';<' and ';>' lines.
  09.09.94 JKB Update to use Read instead of Seq library.
  01.06.07 JKB Supports single-read fasta files; about time too!
  */


#define LINE_LENGTH 60

/*
 * This module should be regarded as part of `read' since it is
 * privy to the internal structure of `Read'.
 *
 * This library also requires use of the mach-io code for the endian
 * independent machine IO.
 *
 * Any references to the writing or reading of edited sequences,
 * or to the bottom strand were added by lfw
 */




/* ---- Imports ---- */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>      /* IMPORT: fopen, fclose, fseek, ftell, fgetc */
#include <ctype.h>      /* IMPORT: isprint */
#include <string.h>

#include "io_lib/misc.h"
#include "io_lib/plain.h"
#include "io_lib/Read.h"
#include "io_lib/xalloc.h"
#include "io_lib/traceType.h"

#include "io_lib/stdio_hack.h"
/* ---- Constants ---- */

#define BasesPerLine 50 /* For output formatting */


/* ---- Exports ---- */


/*
 * Read the plain format sequence from FILE *fp into a Read structure.
 * All printing characters (as defined by ANSII C `isprint')
 * are accepted, but `N's are translated to `-'s.
 *
 * Returns:
 *   Read *     - Success, the Read structure read.
 *   NULLRead   - Failure.
 */
Read *fread_pln(FILE *fp) {
    Read *read = NULLRead;
    off_t fileLen;
    int  ch;
    char *leftc, *rightc, *leftcp, *rightcp;
    int first = 1;

    /*
     * Find the length of the file.
     * Use this as an overestimate of the length of the sequence.
     */
    fseek(fp, (off_t) 0, 2);
    if ((fileLen = ftell(fp)) > INT_MAX /*Was MAXINT2*/)
	goto bail_out;

    fseek(fp, (off_t) 0, 0);
    
    /* Allocate the sequence */
    if (NULLRead == (read = read_allocate(0, fileLen)))
	goto bail_out;

    if (NULL == (leftc = (char *)xmalloc(fileLen)))
	goto bail_out;

    if (NULL == (rightc = (char *)xmalloc(fileLen)))
	goto bail_out;

    leftcp = leftc;
    rightcp = rightc;

    /* Read in the bases */
    
    read->NBases = 0;
    read->format = TT_PLN;

    while ((ch = fgetc(fp)) != EOF) {
	if (ch == '>') {
	    /* Fasta format file - skip the header and load the first
	     * fasta sequence only. We don't even attempt to worry about
	     * multi-sequence file formats for now.
	     */
	    if (!first)
		break;

	    while(ch != '\n' && ch != EOF)
		ch = fgetc(fp);

	}  else if (ch==';') {
	    /*
	     * ;< is left cutoff,
	     * ;> is right cutoff.
	     * Any other ';'s we can treat as a comments.
	     */
	    ch = fgetc(fp);

	    if (first == 1 && ch != '<' && ch != '>') {
		int d;
		char type[5], name[17], line[1024];

		line[0] = ch;
		fgets(&line[1], 1022, fp);

		if (5 == sscanf(line, "%6d%6d%6d%4c%s",
				&d, &d, &d, type, name)) {
		    char * p;

		    if (p = strchr(type, ' '))
			*p = 0;

		    read->format = trace_type_str2int(type);
		    read->trace_name = (char *)xmalloc(strlen(name)+1);
		    if (read->trace_name)
			strcpy(read->trace_name, name);
		}
	    }

	    else if (ch == '<') {
		ch = fgetc(fp);
		while (ch != '\n') {
		    *leftcp++ = ch;
		    ch = fgetc(fp);
		}
	    } else if (ch == '>') {
		ch = fgetc(fp);
		while (ch != '\n') {
		    *rightcp++ = ch;
		    ch = fgetc(fp);
		}
	    } else {
		while(ch != '\n' && ch != EOF)
		    ch = fgetc(fp);
	    }
        } else if (isprint(ch) && !isspace(ch)) {
	    read->base[read->NBases++] = ((ch)=='N') ? '-' : (ch);
	}
	
	first = 0;
    }

    *leftcp = *rightcp = 0;

    read->leftCutoff = strlen(leftc);
    read->rightCutoff = read->leftCutoff + read->NBases + 1;
    memmove(&read->base[read->leftCutoff], read->base, read->NBases);
    memmove(read->base, leftc, read->leftCutoff);
    memmove(&read->base[read->leftCutoff + read->NBases],
	    rightc, strlen(rightc));

    read->NBases += read->leftCutoff + strlen(rightc);
    read->base[read->NBases] = 0;

    xfree(leftc);
    xfree(rightc);
    
    /* SUCCESS */
    return(read);

    /* FAILURE */
 bail_out:
    if (read)
	read_deallocate(read);

    return NULLRead;
}

/*
 * Read the plain format sequence with name `fn' into a Read structure.
 * All printing characters (as defined by ANSII C `isprint')
 * are accepted, but `N's are translated to `-'s.
 *
 * Returns:
 *   Read *     - Success, the Read structure read.
 *   NULLRead   - Failure.
 */
Read *read_pln(char *fn) {
    FILE *fp;
    Read *read;

    /* Open file */
    if ((fp = fopen(fn, "r")) == NULL)
	return NULLRead;

    read = fread_pln(fp);
    fclose(fp);

    if (read && read->trace_name == NULL &&
	(read->trace_name = (char *)xmalloc(strlen(fn)+1)))
	strcpy(read->trace_name, fn);

    return read;
}


/*
 * Write to a Plain file
 */
int fwrite_pln(FILE *fp, Read *read) {
    int i, err = 0;

    for (i = 0; i < read->NBases; i += LINE_LENGTH)
        if (-1 == fprintf(fp, "%.*s\n",
			  read->NBases - i > LINE_LENGTH
			  ? LINE_LENGTH : read->NBases - i,
			  &read->base[i]))
	    err = 1;
    
    return err ? -1 : 0;
}

int write_pln(char *fn, Read *read) {
    FILE *fp;

    if ((fp = fopen(fn,"w")) == NULL) 
	return -1;

    if (fwrite_pln(fp, read)) {
	fclose(fp);
	return -1;
    }

    fclose(fp);
    return 0;
}

