///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BUILD_SCAFFOLD_GRAPH_H
#define BUILD_SCAFFOLD_GRAPH_H

#include "PairsManager.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/Sepdev.h"
#include "paths/reporting/CLinkBundle.h"
#include "paths/reporting/CSuperLinks.h"

/**
 * BuildScaffoldGraph
 *
 * Turn a set of scaffold links into a digraphE. Vertices are oriented
 * scaffolds (0[+], 0[-], 1[+], 1[-], ...), and edges are bundles of
 * links between them, stored as pair<int,int>s of separation and
 * stdev.
 *
 * Filtering on bundles is done in this order:
 *
 * 1. a bundle between two supers is considered at all only if:
 *    . there are at least two links in the bundle, and
 *    . score <= MAX_SCORE, and
 *    . either there are MIN_LINKS links, or the score is >= MIN_SCORE
 *       
 * 2. now look at all the bundles between a super s and another super
 *    s', and between s and the rc of s'. Sort all the bundles by
 *    weight, and look at the weight of the first bundle (w1) against
 *    the weight of the second bundle (w2). If w1 > ( RATIO_MIN_LINKS
 *    - 1 ) * RATIO_TO_SECOND * w2, then we accept the bundle with
 *    weight w1 as the winner of the bundles, and we allow the join
 *    between s and s'fw or s'rc. In this way there is at most one
 *    bundle between any two supers s and s'.
 *
 * 3. we attempt to remove innies with a crude max overlap test: if
 *    the separation between s and s' is < - MAX_OVERLAP and the
 *    weight of the bundle is < LOW_WEIGHT, the bundle is discarded.
 *
 * Other arguments:
 *
 * black_list: if provided, do not allow edges between these supers
 * slop: argument to AllLinks( ) in class CSuperLink
 */
void BuildScaffoldGraph( const PairsManager &pairs,
			 const vec<superb> &supers,
			 const vec<alignlet> &aligns,
			 const vec<int> &index,
			 digraphE<sepdev> &graph,
			 digraphE<CLinkBundle> *bgraph = 0,
			 vec< pair<int,int> > *black_list = 0,
			 ostream *log = 0,
			 float MAX_SCORE = 1.65,
			 int MIN_LINKS = 4,
			 float MIN_SCORE = 0.5,
			 int RATIO_MIN_LINKS = 2,
			 double RATIO_TO_SECOND = 6.0,
			 int MAX_OVERLAP = 10000,
			 int LOW_WEIGHT = 6,
			 double slop = 3.5 );

#endif
