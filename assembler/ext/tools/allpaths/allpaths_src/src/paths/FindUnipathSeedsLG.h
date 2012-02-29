///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef FIND_UNIPATH_SEEDS_LG_H
#define FIND_UNIPATH_SEEDS_LG_H

// FindUnipathSeeds: pick seed unipaths for LocalizeReads.

#include "CoreTools.h"
#include "PairsManager.h"
#include "ReadLocationLG.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerPath.h"
#include "paths/simulation/Placement.h" // placement
#include "paths/UnipathNhoodLG.h" // sepdev, fsepdev



enum SeedStatus
  { SEED_GOOD, SEED_ISOLATED, SEED_SHORT, SEED_HIGH_CN, SEED_ZERO_CN, SEED_RC_ON_REF, SEED_REDUNDANT, SEED_RC_OF_SEED };

static const string status_names[] =
  {"=== GOOD ===", "isolated", "short", "high-copy-number", "zero-copy-number", "rc-on-reference", "redundant", "rc-of-another-seed"};



void FindUnipathSeeds( 

     // inputs:

     const int K,
     const int ploidy,
     const int MIN_KMERS_IN_SEED,      // smallest seed to consider
     const int MIN_KMERS_IN_NEIGHBORS, // seed plus neighbors must be this big
     const vec<int>& ulen,             // unipath lengths
     const vec<unipath_id_t>& to_rc,   // unipath involutions
     const vec<int>& predicted_copyno, // predicted copy number for unipaths
     const digraphE<fsepdev>& FG,      // graph of all normal unipaths
     const Bool USE_TRUTH,             // if so, cheat by discarding RC unipaths
     VecPlacementVec locs,             // unipath locations on reference
     const int MAX_SEED_DIST,          // must use seed if farther than this
     const vec<Bool> & branches, 
     const vecKmerPath* global_paths,
     const vecKmerPath* global_paths_rc,
     const vec<tagged_rpint>* global_pathsdb,
     const vec<int>* global_read_lengths,
     const vec<ReadLocationLG>* global_unilocs,
     const vec<longlong>* global_unilocs_index, 
     const PairsManager* global_pairs,

     // output:

     vec<int>& seeds,
     vec<SeedStatus>& seed_status,
     ostream & seed_log,

     // extra params:
     
     const unsigned NUM_THREADS
          );



void
EvalUnipathSeeds( const vec<int> & seeds, const vec<int>& ulen,
		  const vecbasevector & genome, const VecPlacementVec& locs,
		  const int MAX_SEED_DIST, const vec<int> & predicted_copyno,
		  const digraphE<fsepdev> & FG, const vec<SeedStatus>& seed_status );

#endif
