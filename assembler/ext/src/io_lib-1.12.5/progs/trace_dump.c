/*
 * Copyright (c) Medical Research Council 2002. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written by James Bonfield, Simon Dear, Rodger Staden,
 * Mark Jordan as part of the Staden Package at the MRC Laboratory of
 * Molecular Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <io_lib/Read.h>

int main(int argc, char **argv)
{
    Read* read;
    int i;

    if (argc != 2) {
    fprintf(stderr, "Usage: trace_dump <trace file>\n");
    return 1;
    }


    read = read_reading( argv[1], TT_ANY );


    if (read == NULL) {
    fprintf(stderr, "Tracedump was unable to open file %s\n", argv[1] );
    return 1;
    }

    printf("[Trace]\n");
    printf("%s\n", read->trace_name );

    printf("\n[Header]\n");
    printf("%d\t\t# format\n",          read->format);
    printf("%d\t\t# NPoints\n",         read->NPoints);
    printf("%d\t\t# NBases\n",          read->NBases);
    printf("%d\t\t# NFlows\n",          read->nflows);
    printf("%d\t\t# maxTraceVal\n",     (int)read->maxTraceVal-read->baseline);
    printf("%d\t\t# baseline\n",        read->baseline);
    printf("%d\t\t# leftCutoff\n",      read->leftCutoff);
    printf("%d\t\t# rightCutoff\n",     read->rightCutoff);

    puts("\n[Bases]");
    for (i = 0; i < read->NBases; i++) {
    printf("%c %05d %+03d %+03d %+03d %+03d #%3d\n",
           read->base[i],
           read->basePos ? read->basePos[i] : 0,
           (int)read->prob_A[i],
           (int)read->prob_C[i],
           (int)read->prob_G[i],
           (int)read->prob_T[i],
           i);
    }

    if (read->NPoints) {
	puts("\n[A_Trace]");
	for(i = 0; i < read->NPoints; i++)
	    printf("%d\t#%5d\n", (int)read->traceA[i] - read->baseline, i);

	puts("\n[C_Trace]");
	for(i = 0; i < read->NPoints; i++)
	    printf("%d\t#%5d\n", (int)read->traceC[i] - read->baseline, i);

	puts("\n[G_Trace]");
	for(i = 0; i < read->NPoints; i++)
	    printf("%d\t#%5d\n", (int)read->traceG[i] - read->baseline, i);

	puts("\n[T_Trace]");
	for(i = 0; i < read->NPoints; i++)
	    printf("%d\t#%5d\n", (int)read->traceT[i] - read->baseline, i);
    }

    if (read->flow_order) {
	puts("\n[Flows]");
	for (i = 0; i < read->nflows; i++) {
	    printf("%c %5.2f  %u\t#%5d\n",
		   read->flow_order[i],
		   read->flow ? read->flow[i] : 0,
		   read->flow_raw ? read->flow_raw[i] : 0,
		   i);
	}
    }

    if (read->info) {
	puts("\n[Info]");
	printf("%s\n", read->info);
    }

    read_deallocate(read);

    return 0;
}
