 /////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#ifndef MAKE_SCAFFOLDS_CLOSE_BEST_H
#define MAKE_SCAFFOLDS_CLOSE_BEST_H

#include "Fastavector.h"
#include "PairsManager.h"
#include "Superb.h"
#include "VecUtilities.h"
#include "paths/Alignlet.h"

/**
 * MakeScaffoldsCloseBest
 *
 * Scaffold contigs together using aligments of reads onto the
 * contigs. It first generates a scaffold graph (see below), and then
 * it iteratively joins oriented scaffolds by selecting the closest
 * best next scaffold.
 *
 * The graph's vertices are, in this order: s0fw, s0rc, s1fw,
 * s1rc,... 
 *
 * MIN_PAIR_SEPS
 * MAX_PAIR_SEPS
 * MIN_LINKS: these three arguments must be in sync. For each
 *   iteration, accept links between two scaffolds if and only if
 *   there are at least MIN_LINKS links of size in the interval
 *   specified by MIN_PAIR_SEPS and MAX_PAIR_SEPS. Notice: if
 *   MIN_PAIR_SEPS is empty, both MIN_PAIR_SEPS and MAX_PAIR_SEPSD
 *   will be derived from library information.
 * MAX_LINKS: similar to MAX_LINKS (if not empty). If given, it must
 *   be in sycn with the three args above.
 *
 * MAX_OVERLAP: exclude links implying a large overlap
 * MIN_SCAFFOLD_LEN: only allow links between long scaffolds
 * 
 * SUCK_SCAFFOLDS: if true, run code to absorb small scaffolds
 * SHAVE_GRAPH: if true, shave the scaffold graph
 *
 * SCAFFOLD_GRAPH_OUT: if not empty, save scaffold graphs (one per iteration)
 * SCAFFOLD_GRAPH_DOT: if not empty, save scaffold graphs as dot file (one
 *   per iteration)
 *
 * VERBOSITY: verbosity flag (0 to 3, 0 is less verbose, 3 is a deluge)
 * MIN_LINKS_TO_PRINT: logging args
 */
void MakeScaffoldsCloseBest ( vec<superb> &scaffolds,
			      vec<fastavector> &scaffold_contigs,
			      vec<alignlet> &aligns0,
			      vec<int> &aligns0_index,
			      const PairsManager &pairs,
			      
			      const String MIN_PAIR_SEPS = "",
			      const String MAX_PAIR_SEPS = "",
			      const String MIN_LINKS = "{20,6,4,2}",
			      const String MAX_LINKS = "",

			      const int MAX_OVERLAP = std::numeric_limits<int>::max(),
			      const int MIN_SCAFFOLD_LEN = 0,
			      
			      const Bool SUCK_SCAFFOLDS = False,
			      const Bool SHAVE_GRAPH = False,
			      
			      const String SCAFFOLD_GRAPH_OUT = "",
			      const String SCAFFOLD_GRAPH_DOT = "",
			      
			      const int VERBOSITY = 0,
			      const int MIN_LINKS_TO_PRINT = std::numeric_limits<int>::max( ) );

#endif
