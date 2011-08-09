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
 * File: 	read_alloc.c
 * Purpose:	Performs the allocation/freeing of Read structures
 * Last update: 01/09/94
 */


/*
    The Read data type is designed so that it can hold a varying degree
    of information about sequences, yet have a single set of calls
    to access the data.

    There are plenty of assumptions around that both the number of
    bases and the number of points will fit into an int_2, a short.

*/

/* ---- Includes ---- */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "io_lib/misc.h"
#include "io_lib/Read.h"

/*
 * Allocate a new sequence, with the given sizes.
 * Returns:
 *   "Read *" for success
 *   "NULLRead" for failure
 */
Read *read_allocate(int num_points, int num_bases) {
    Read *seq = NULLRead;

    int sections = read_sections(0);

    /* Allocate the body of the sequence */
    if ((seq = (Read *)xmalloc(sizeof(Read))) == NULL)
	return(NULLRead);

    seq->NPoints = num_points;
    seq->NBases  = num_bases;

    /*   
     * Initialise the body, all pointers are set to NULL so we can
     * happily call `read_deallocate()`.
     */
    seq->leftCutoff  = 0;
    seq->rightCutoff = 0;
    seq->maxTraceVal = 0;
    seq->baseline = 0;

    seq->traceC    = NULL;
    seq->traceA    = NULL;
    seq->traceG    = NULL;
    seq->traceT    = NULL;

    seq->base      = NULL;
    seq->basePos   = NULL;

    seq->info = NULL;
    seq->format = TT_ANY;
    seq->trace_name = NULL;

    seq->prob_A = NULL;
    seq->prob_C = NULL;
    seq->prob_G = NULL;
    seq->prob_T = NULL;

    seq->orig_trace_format = TT_ANY;
    seq->orig_trace = NULL;
    seq->orig_trace_free = NULL;

    seq->ident = NULL;

    /* Allocate space for the bases - 1 extra for the ->base field so
     * that we can treat it as a NULL terminated string.
     */
    if (sections & READ_BASES &&
	(((seq->base	  = (char *)xcalloc(num_bases+1,1))   == NULL) ||
	 ((seq->basePos   = (uint_2 *)xcalloc(num_bases+1,2)) == NULL) ||
	 ((seq->prob_A    = (char *)xcalloc(num_bases+1,1))   == NULL) ||
	 ((seq->prob_C    = (char *)xcalloc(num_bases+1,1))   == NULL) ||
	 ((seq->prob_G    = (char *)xcalloc(num_bases+1,1))   == NULL) ||
	 ((seq->prob_T    = (char *)xcalloc(num_bases+1,1))   == NULL))
	)
    {
	read_deallocate(seq);
	return NULLRead;
    }

    if (sections & READ_SAMPLES &&
	(((seq->traceC   =(TRACE *)xcalloc(num_points+1, 2))  == NULL)||
	 ((seq->traceA   =(TRACE *)xcalloc(num_points+1, 2))  == NULL)||
	 ((seq->traceG   =(TRACE *)xcalloc(num_points+1, 2))  == NULL)||
	 ((seq->traceT   =(TRACE *)xcalloc(num_points+1, 2))  == NULL))
	)
    {
	read_deallocate(seq);
	return NULLRead;
    }

    seq->nflows = 0;
    seq->flow_order = NULL;
    seq->flow = NULL;
    seq->flow_raw = NULL;

    seq->private_data = NULL;
    seq->private_size = 0;
    
    return seq;
}


/*
 * Free memory allocated to a sequence by read_allocate().
 */
void read_deallocate(Read *read)
{
    if (read == NULLRead)
	return;

    if (read->traceC  != NULL)  xfree(read->traceC);
    if (read->traceA  != NULL)  xfree(read->traceA);
    if (read->traceG  != NULL)  xfree(read->traceG);
    if (read->traceT  != NULL)  xfree(read->traceT);

    if (read->base    != NULL)  xfree(read->base);
    if (read->basePos != NULL)  xfree(read->basePos);

    if (read->info    != NULL)  xfree(read->info);

    if (read->prob_A  != NULL)  xfree(read->prob_A);
    if (read->prob_C  != NULL)  xfree(read->prob_C);
    if (read->prob_G  != NULL)  xfree(read->prob_G);
    if (read->prob_T  != NULL)  xfree(read->prob_T);

    if (read->trace_name != NULL) xfree(read->trace_name);

    if (read->orig_trace != NULL) {
	if (read->orig_trace_free)
	    read->orig_trace_free(read->orig_trace);
	else
	    xfree(read->orig_trace);
    }

    if (read->ident != NULL)
	xfree(read->ident);

    if (read->flow_order)
	xfree(read->flow_order);
    if (read->flow)
	xfree(read->flow);
    if (read->flow_raw)
	xfree(read->flow_raw);

    if (read->private_data)
	xfree(read->private_data);

    xfree(read);
}




/*
 * Duplicates the read structure and optionally gives it a new filename.
 * The following fields are not duplicated:
 *    
 *  int  orig_trace_format;
 *  void (*orig_trace_free)(void *ptr);
 *  void *orig_trace;
 *  char *ident;
 *
 * Returns:
 *   "Read *" for success
 *   "NULLRead" for failure
 */
Read* read_dup( Read* src, const char* new_name )
{
    int   n;
    Read* dst;
    assert(src);

    /* Allocate storage and initialise */
    dst = read_allocate( src->NPoints, src->NBases );
    if( dst == NULLRead )
	return 0;
    dst->info       = 0;
    dst->trace_name = 0;


    /* Copy over possibly new name */
    if( new_name )
	n = strlen(new_name);
    else if( src->trace_name )
	n = strlen(src->trace_name);
    else
	n = 0;
    if( n > 0 )	{
	dst->trace_name = (char*) xmalloc(n+1);
	if( !dst->trace_name )
	    goto error;

	if(new_name) 
	    strcpy( dst->trace_name, new_name );
	else
	    strcpy( dst->trace_name, src->trace_name );
    }
	
	
    /* Copy over info */
    if( src->info ) {
	dst->info = strdup(src->info);
    }


    /* Copy single fields */
    dst->format      = src->format;
    dst->maxTraceVal = src->maxTraceVal;
    dst->leftCutoff  = src->leftCutoff;
    dst->rightCutoff = src->rightCutoff;
    dst->baseline    = src->baseline; 
    

    /* Copy NPoints fields if they exist */
    if( src->traceA )
	{
	    for( n=0; n<src->NPoints; n++ )
		{
		    dst->traceA[n] = src->traceA[n];
		    dst->traceC[n] = src->traceC[n];
		    dst->traceG[n] = src->traceG[n];
		    dst->traceT[n] = src->traceT[n];
		}
	}
    
    
    /* Copy NBases fields if they exist */
    if( src->base && src->base[0] )
	{
	    for( n=0; n<src->NBases; n++ )
		{
		    dst->base[n]    = src->base[n];
		    dst->basePos[n] = src->basePos[n];
		    if( src->prob_A )
			{
			    dst->prob_A[n] = src->prob_A[n];
			    dst->prob_C[n] = src->prob_C[n];
			    dst->prob_G[n] = src->prob_G[n];
			    dst->prob_T[n] = src->prob_T[n];
			}
		}
	}
    
    
    /* Success */
    return dst;

 error:
    /* Failure */
    read_deallocate(dst);
    return NULLRead;
}
