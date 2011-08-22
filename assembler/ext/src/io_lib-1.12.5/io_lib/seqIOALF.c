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
 * Title:       seqIOALF
 *
 * File: 	 seqIOALF.c
 * Purpose:	 IO of ALF sequences
 * Last update:  9th September 1994
 */

/*
 * Change Log :- 
 * 14.01.91 SD
 * when complimenting the sequence with an odd number of bases,
 * the middle base position was not adjusted.
 * 15.01.91 SD  Put StLouis stuff on compilation flag
 * 15.01.91 SD  New include file (opp.h)
 * 02.08.91 SD  Changes the mapping of uncertainty codes so that we
 * now only generate A C G T and -
 * Previously... bug in interpreting ALF integer fields.
 * We now treat them as unsigned.
 * 17.09.91 LFW changed STLOUIS compilation flag to SAVE_EDITS
 * and AUTO_CLIP
 * 25.10.91 SD  Machine independant I/O...removed BIGENDIAN flag
 * 25.11.91 SD There was a hard limit (of 1024) for allocation of
 * space for number of bases, yet program would 
 * read in more if there were any, causing nasties to happen.
 * 
 * 11.11.92 LFW added section to actually check that the trace it
 * is trying to open is an ALF file using traceType sub
 * 
 * 10.11.92 SD  SCF comments now stored in seq data structure
 * 09.09.94 JKB Update to use Read instead of Seq library.
 * 04.03.98 JKB Look for "Raw data" when "Processed data" is not found.
 */

/* RMD I made substantial changes to this file 12/28/90 so as to
 * read sequence data more freely (necessary when reading data from
 * multiple trace files).
 * The affected area is indicated by comments starting RMD, like
 * this one.
 */

/* This file was adapted by LFW from seqIOABI.c.
 * The ALF results file is a concatenation of many files with an
 * index structure at the beginning, consisting of a 512 byte
 * block that we ignore, followed by 128 byte blocks describing
 * each file.  All files, including the header region, are rounded 
 * up to a multiple of 512 bytes long.  
 * The getIndexEntry routines identify the 128 byte index component
 * of interest by matching 4 chars of its ASCII label, then extract
 * the field of choice from that entry.
 * 
 * Note that the SUN and PC are of opposite endian-ness, so that
 * we have to provide special routines to read words and longwords
 * from the results file.  Luckily the floating point numbers are
 * written out in ASCII.
 */


/* ---- Imports ---- */


#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "io_lib/stdio_hack.h"

#include "io_lib/Read.h"
#include "io_lib/mach-io.h"
#include "io_lib/xalloc.h"

/* ---- Constants ---- */

#define BasesPerLine 50 /* For output formatting */

#define IndexEntryLength ((off_t)128)


/*
 * Here are some labels we will be looking for, four chars packed
 * into a long word.
 */
#define EntryLabel        ((uint_4) ((((('A'<<8)+'L')<<8)+'F')<<8)+' ')
#define BaseEntryLabel    ((uint_4) ((((('S'<<8)+'e')<<8)+'q')<<8)+'u')
#define DataEntryLabel    ((uint_4) ((((('P'<<8)+'r')<<8)+'o')<<8)+'c')
#define RawDataEntryLabel ((uint_4) ((((('R'<<8)+'a')<<8)+'w')<<8)+' ')

/* RMD make enough space for bases - hard limit */
#define BASELIMIT 4096


/* ---- Internal functions ---- */

/*
 * From the ALF results file connected to `fp' whose index starts
 * at byte offset `indexO', return in `val' the `lw'th long word
 * from the entry labelled `label'.
 * The result is 0 for failure, 1 for success.
 */
static int getIndexEntryLW(FILE *fp, off_t indexO,
			   uint_4 label, int lw,
			   uint_4 *val) {
    off_t entryNum=-1;
    int i;
    uint_4 entryLabel;
    
    do {
	entryNum++;
	if (fseek(fp, indexO+(entryNum*IndexEntryLength), 0) != 0) 
	    return 0;
	    
	if (!be_read_int_4(fp, &entryLabel))
	    return 0;
    } while (!(entryLabel == label));
    
    for(i=2; i<lw; i++)
	if (!be_read_int_4(fp, val))
	    return 0;
    
    
    /* when i = lw read in the 4 bytes backwards */
    if (!le_read_int_4(fp,val))
	return 0;
    
    return 1;
}

/*
 * From the ALF results file connected to `fp' whose index starts
 * at byte offset `indexO', return in `val' the `lw'th  word (int2)
 * from the entry labelled `label'.
 * The result is 0 for failure, 1 for success.
 */
static int getIndexEntryW(FILE *fp, off_t indexO,
			  uint_4 label, int lw,
			  uint_2 *val) {
    off_t entryNum=-1;
    int i;
    uint_4 entryLabel;
    uint_4 jval;
    
    do {
	entryNum++;
	if (fseek(fp, indexO+(entryNum*IndexEntryLength), 0) != 0)
	    return 0;

	if (!be_read_int_4(fp, &entryLabel))
	    return 0;
	} while (!(entryLabel == label));
    
    
    for(i=2; i<lw; i++)
	if (!be_read_int_4(fp, &jval))
	    return 0;

    if (!le_read_int_2(fp, val))
	return 0;
    
    return 1;
}


/* ---- Exports ---- */


/*
 * Read the ALF format sequence from FILE *fp into a Read structure.
 * All printing characters (as defined by ANSII C `isprint')
 * are accepted, but `N's are translated to `-'s. In this respect we
 * are adhering (more or less) to the CSET_DEFAULT uncertainty code set.
 * 
 * Returns:
 *   Read *	- Success, the Read structure read.
 *   NULLRead	- Failure.
 */
Read *fread_alf(FILE *fp) {
    Read *read = NULLRead;
    int i;
    int numPoints;
    int sections = read_sections(0);
    
    uint_4 data_size;
    uint_4 dataO;
    uint_4 header_size=396; /* size of the header of the processed data
			       section */
    uint_2 actBaseDataSize; /* actual number of bytes of data of information
			       containing the base and basePos information */
    int num_points;         /* keeps track of the actual number of points,
			       rather than the early guess of numPoints */

    off_t indexO;           /* File offset where the index is */
    uint_4 baseO;           /* File offset where the bases are stored */
    
    
    /*
     * RMD lots of changes below here until end of data reading section
     * Some are cosmetic.
     * getIndexEntry calls in front of where they were needed, and made
     * There is a substantive change to the inner loop of the sequence
     * reading section.  This now uses fscanf - much less rigid than the
     * previous scheme.  Note that it reads bp as a float.  This is because
     * it is a float in multiple trace data files! (bizarre Pharmacia
     * programming!).
     */
    
    
    /*************************************************************
     * Read the various file offsets
     *************************************************************/

    /* indexO is the offset of the index.
     * Or I could look for the first label, starting 'ALF'
     * if I used 512 then none of the entries are on long 
     * word boundaries
     */
    indexO = 522;
    
    /* offset in file of first base of sequence */
    if (! (getIndexEntryLW(fp,indexO,BaseEntryLabel,12,&baseO)) )
	goto bail_out;
    
    /* actual size of region containing this data */
    if (! (getIndexEntryW(fp,indexO,BaseEntryLabel,10,&actBaseDataSize)) )
	goto bail_out;
    
    /* Look for Processed data first. If we fail to find it, then look for
     * the Raw data (same format).
     */

    /* offset in file to start of processed data segment - there 
     * is then a header of size header_size (currently 396)
     */
    if (! (getIndexEntryLW(fp,indexO,DataEntryLabel,12,&dataO)) ) {
	if (! (getIndexEntryLW(fp,indexO,RawDataEntryLabel,12,&dataO)) )
	    goto bail_out;

	/* actual size of region containing this data */
	if (! (getIndexEntryLW(fp,indexO,RawDataEntryLabel,10,&data_size)) )
	    goto bail_out;
    } else {
	/* actual size of region containing this data */
	if (! (getIndexEntryLW(fp,indexO,DataEntryLabel,10,&data_size)) )
	    goto bail_out;
    }
    
    /* Because each trace value is stored in a 2 byte
     * integer, thus to store A C G T information
     * it takes 8 bytes.  So subtract off the header and
     * divide by 8
     */
    numPoints = (int)((data_size - header_size)/ 8); 
    
    /* Allocate the sequence */
    if (NULLRead == (read = read_allocate(numPoints, BASELIMIT)))
	goto bail_out;
    
    /*************************************************************
     * Read the bases information
     *************************************************************/
    if (sections & READ_BASES) {
	/* new locals introduced by LFW and/or RMD for the ALF */
	int numBases;	/* number of nucleotides read in */
	float bp;
	char ch;

	if (!(fseek(fp, (off_t)baseO, 0) == 0))
	    goto bail_out;
	
	for (numBases = 0; (unsigned)ftell(fp) < baseO+(unsigned short)actBaseDataSize
	                   && numBases<BASELIMIT;) {
	    char line[200];

	    fgets(line, (int)sizeof(line), fp);
	    sscanf(line, "%c %*d %f", &ch, &bp);
	    
	    /* we convert ch to Staden format here */
	    switch (ch) {
	    case 'A':
	    case 'C':
	    case 'G':
	    case 'T':
		break;
	    default:
		ch = '-';
/*
		if (isupper(ch))
		    ch = '-';
		else
		    ch = '\0';
*/
	    }
	    
	    if (ch) {
		read->base[numBases]    = ch;
		read->prob_A[numBases]	= 0;
		read->prob_C[numBases]	= 0;
		read->prob_G[numBases]	= 0;
		read->prob_T[numBases]	= 0;
		read->basePos[numBases] = bp;
		++numBases;
	    }
	}
	read->base[numBases] = 0;
	
	read->NBases  = numBases;
    }
    
    /*************************************************************
     * Read the trace information
     *************************************************************/
    
    if (sections & READ_SAMPLES) {
	
	/*
	 * Traces are stored as 2 byte integers in records in the order of
	 * A C G T A C G T ...
	 */
	
	if (fseek(fp, (off_t)(dataO+header_size), 0) != 0) 
	    goto bail_out;
	
	num_points = 0;
	
	for (i=0; i < read->NPoints; i++) {
	    if (!le_read_int_2(fp, &(read->traceA[i])))
		goto bail_out;
	    if (read->maxTraceVal < read->traceA[i])
		read->maxTraceVal = read->traceA[i];
	    
	    if (!le_read_int_2(fp, &(read->traceC[i])))
		goto bail_out;
	    if (read->maxTraceVal < read->traceC[i])
		read->maxTraceVal = read->traceC[i];
	    
	    if (!le_read_int_2(fp, &(read->traceG[i])))
		goto bail_out;
	    if (read->maxTraceVal < read->traceG[i])
		read->maxTraceVal = read->traceG[i];
	    
	    if (!le_read_int_2(fp, &(read->traceT[i])))
		goto bail_out;
	    if (read->maxTraceVal < read->traceT[i])
		read->maxTraceVal = read->traceT[i];
	    
	    if (read->traceA[i]==0 && read->traceT[i]==0 &&
		read->traceC[i]==0 && read->traceG[i]==0 &&
		i > (numPoints-64))
		break;
	    
	    num_points++;
	}
    }
    
    /* SUCCESS */

    read->format = TT_ALF;
    return(read);

    /* FAILURE */
 bail_out:
    if (read)
	read_deallocate(read);

    return NULLRead;
}

/*
 * Read the ALF format sequence with name `fn' into a Read structure.
 * All printing characters (as defined by ANSII C `isprint')
 * are accepted, but `N's are translated to `-'s. In this respect we
 * are adhering (more or less) to the CSET_DEFAULT uncertainty code set.
 * 
 * Returns:
 *   Read *	- Success, the Read structure read.
 *   NULLRead	- Failure.
 */
Read *read_alf(char *fn) {
    FILE *fp;
    Read *read;

    /* Open file */
    if ((fp = fopen(fn, "rb")) == NULL)
	return NULLRead;

    read = fread_alf(fp);
    fclose(fp);

    if (read && (read->trace_name = (char *)xmalloc(strlen(fn)+1)))
	strcpy(read->trace_name, fn);

    return read;
}

/*
 * Write to an ALF file - unsupported.
 */
/* ARGSUSED */
int write_alf(char *fn, Read *read) {
    fprintf(stderr, "ALF write support is unavailable\n");
    return -1;
}

/*
 * Write to an ALF file - unsupported.
 */
/* ARGSUSED */
int fwrite_alf(FILE *fp, Read *read) {
    fprintf(stderr, "ALF write support is unavailable\n");
    return -1;
}
