///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__FLATTEN_HYPER_FASTAVECTOR_H
#define PATHS__FLATTEN_HYPER_FASTAVECTOR_H

#include "paths/CLinFastavec.h"
#include "paths/SquashHyperFastavector.h"

/**
 * FlattenHyperFastavector
 *
 * hfv: will be modified
 * base_out: save output only if not empty!
 * hyper_inter_file: if not empty, save intermediate hypers
 * dump_hfv_head: if not empty, save hyper before and after flattening
 * INITIAL_SCAFFOLDS_PER_CONTIG: break scaffolds into singleton scaffolds
 * NUM_THREADS: for FirstLookupFinder
 */
void FlattenHyperFastavector( ostream &log,
			      HyperFastavector &hfv,
			      const String base_out = "",
			      const String hyper_inter_file = "",
			      const String dump_hfv_head = "", 
			      const bool NEW_ALGORITHM = False,
			      const bool INITIAL_SCAFFOLD_PER_CONTIG = True,
			      const int MAX_CELL_SIZE = 20,
			      const int MIN_EDGE_TO_SAVE = 1000,
			      const int NUM_THREADS = 0 );

#endif
