///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef UNIPATH_NHOOD_LG_H
#define UNIPATH_NHOOD_LG_H

// File: UnipathNhoodLG.h
//
// This file defines a toolkit for building a neighborhood (abbreviated "nhood")
// of unipaths around a given seed unipath, and identifying the reads that came from
// this neighborhood.

#include "CoreTools.h"
#include "Equiv.h"
#include "Basevector.h"
#include "Intvector.h"
#include "PairsManager.h"
#include "ReadLocationLG.h"
#include "SemanticTypes.h"
#include "VecUtilities.h"
#include "graph/Digraph.h"
#include "paths/KmerPath.h"
#include "paths/PdfEntry.h"
#include "paths/simulation/Placement.h"
#include "paths/UnipathNhoodLG.h"
#include "paths/UnipathNhoodCommon.h"

// FuncDecl: BuildUnipathLinkGraph
//
// Build the unipath graph, in which vertices are 
// normal unipaths and edges come from read pairs.
//
// Instantiated for sepdev and fsepdev in the .cc file.
template<class T>
void BuildUnipathLinkGraph( 

     // inputs:

     const int K,                              // as in Kmer
     const int ploidy,
     const PairsManager& pairs,                // read pairs
     const vec<ReadLocationLG>& ulocs,         // locations of reads on unipaths
     const VecULongVec& ulocs_indexr,          // index to it by reads
     const vec<Bool>& normal,                  // is a given unipath normal
     const vec<int>& ulen,                     // unipath lengths
     const vec<int>& rlen,                     // read lengths (in kmers)
     const vec<unipath_id_t>& to_rc,           // map unipath to its rc
     const int min_edge_multiplicity,          // else ignore edge
     const vec<int>& CNs,                      // unipath copy numbers

     // output:

     digraphE< Tsepdev<T> >& G,                // the graph

     // verbosity:

     Bool verbose = False );

/// FuncDecl: FillInTransitiveEdges
///
/// For each vertex of predicted-copy-number one in the graph, join it
/// to other vertices whose separation will be within radius of that
/// vertex.  But this naive version won't find connections which
/// require passing though multiple vertices further than distance
/// radius away.
///
/// An edge we construct will replace an existing edge with the same
/// endpoints if the new one has smaller deviation.
///
/// Do not introduce new edges having deviation > max_dev.
///
/// Do not replace an edge unless its deviation decreases by percent_improvement.
template<class T>   // defined for int and double
void FillInTransitiveEdges( digraphE< Tsepdev<T> >& graph, 
			    const int radius, 
                            const double max_dev,
                            const double percent_improvement,
			    const vec<int>& predicted_copyno,
                            const int ploidy,
			    const vec<nbases_t>& unipath_len, 
                            const int min_unipath = 0 );













#endif
