///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef CORRECT_LONG_READS_1_H
#define CORRECT_LONG_READS_1_H

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/CorrectLongReadsTools.h"
#include "paths/LongReadTools.h"
#include "paths/Uniseq.h"

void Phase1( 

     // inputs:

     const int NUM_THREADS, const vecbasevector& unibases, const vec<int>& to_rc,
     const vec< vec<int> >& nexts, const int K, const int L,
     const vec< vec< pair<int,int> > >& Ulocs, const vecbasevector& longreads,
     const vec<int>& reads_to_use, const heuristics& heur, 

     // algorithms:

     const Bool USE_SHORTEST, const Bool CLUSTER_ALIGNS_NEW_CLEAN,

     // control:

     const Bool PATCHES_ONLY, const int MIN_PATCH1, const int MIN_PATCH2,
     const Bool STANDARD_ALIGNS, const Bool FILTER,
     const Bool CLEAN_GRAPH, const Bool CYCLIC_BAIL,
     const Bool KILL_INFERIOR_CLUSTERS_NEW, const int MIN_SPREAD1,

     // logging:

     const vec<int>& reads_to_trace, const Bool VERBOSE1, const Bool SKIP_SILENT, 
     const Bool DOT1, const Bool QLT1, const String& data_dir, 
     const String& run_dir, const Bool ABBREVIATE_ALIGNMENTS,
     const Bool PRINT_RAW_ALIGNS,

     // outputs:

     vec<uniseq>& UNISEQ, vec< vec<int> >& UNISEQ_ID, vec<Bool>& COMPUTED, 
     vecbasevector& all, vec<GapPatcher>& patchers, 
     vec< digraphVE<int,int> >& Hall );

#endif
