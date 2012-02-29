/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef MAKE_SCAFFOLDS_SCORED_BEST_H
#define MAKE_SCAFFOLDS_SCORED_BEST_H

#include "Fastavector.h"
#include "PairsManager.h"
#include "Superb.h"
#include "Vec.h"
#include "paths/Alignlet.h"

/**
 * MakeScaffoldsScoredBest
 *
 * Make scaffolds by generating a scaffold graph (see below) and
 * iteratively joining oriented supers by selecting first edges with
 * the best edge score.
 *
 * The graph's vertices are, in this order: s0fw, s0rc, s1fw,
 * s1rc,... The graphs edges are CLinkBundle objects, ie bundles of
 * consistent links between oriented supers. The edges are sorted by
 * score, and the best edges are used to iteratively join vertices.
 *
 * The following tests are run only if the arg pointers are not NULL,
 * in this order:
 *   MIN_LINKS: if given, discard edges with a weight < *MIN_LINKS
 *   MIN_SEP: if given, discard edges with separation < *MIN_SEP
 *
 * SUCK_SCAFFOLDS: run SuckScaffolds (aligns are flipped!)
 * VERBOSE: log verbosity
 */
void MakeScaffoldsScoredBest( const PairsManager &pairs,
			      vec<fastavector> &contigs,
			      vec<superb> &supers,
			      vec<alignlet> &aligns,
			      vec<int> &index,
			      ostream &out,
			      int *MIN_LINKS = 0,
			      int *MIN_SEP = 0,
			      bool SUCK_SCAFFOLDS = true,
			      bool VERBOSE = false );

#endif
