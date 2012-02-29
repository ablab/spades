/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// This code takes a set of reads and calculates a set of "super
// reads" by finding the unique sequences of unipaths into which they
// fall.  Finding the muxes among this set of super reads is generally
// much faster than among all the reads.

// Optionally, we also calculate all the possible distances of each
// unique unipath sequence from the others in the resulting mux graph;
// this will allow us to identify whether some mux is a valid distance
// away from the target read when we're searching.

// We add in trivial muxes to tell us where the given reads fall in
// the unique unipaths that contain them, and generate the subsumption
// list data as well.  We convert the unique unipath sequences back in
// the KmerPaths and append those to the original read set.

#include "ReadPairing.h"

#include "paths/KmerPath.h"
#include "paths/MuxGraph.h"
#include "paths/OffsetTracker.h"
#include "paths/SubsumptionList.h"

void AddSuperReads( const vecKmerPath& pathsFw, 
                    const vecKmerPath& pathsRc, 
                    const vec<read_pairing>& pairs,
                    const int K, 
                    const int MAXSEP, 
                    const Float sdMult,
                    vecKmerPath& allPathsFw, 
                    vecKmerPath& allPathsRc,
                    MuxGraph& allMuxes, 
                    SubsumptionList& allSubs,
                    OffsetTracker* pTracker );
