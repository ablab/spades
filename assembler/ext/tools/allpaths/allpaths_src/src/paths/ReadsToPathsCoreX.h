/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/* ReadsToPathsCoreX: The pathing function.
 *
 * Take an input vecbasevector and path it with the given value of K, producing
 * a vecKmerPath.
 *
 *
 * Parallelization is controlled via the PARCEL_DIR and NUM_THREADS options.
 * A set of KmerParcels will be created with parallelization level NUM_THREADS,
 * which should not be more than the number of CPUs this process is to use.
 *
 * If a PARCEL_DIR is specified, the KmerParcels will be stored in the given
 * directory via a KmerParcelsDiskStore; otherwise they will be stored in local
 * memory.  Storing on disk reduces memory overhead, but it increases I/O use;
 * it is also unsafe in the event of multiple simultaneous calls to
 * ReadsToPathsCoreY from parallel threads.
 *
 * Use of the CHECKPOINT option allows you to stop a pathing process and pick
 * up where you left off later; see MakeAlignsPathsParallelX.cc for details.
 *
 *****************************************************************************/

#ifndef READS_TO_PATHS_COREX_H
#define READS_TO_PATHS_COREX_H

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/KmerPath.h"


// NOTE: 'bvv' is not 'const' because we need to sort it.
//   At the end, 'bvv' is returned to it's original ordering.
//   The alternative is to keep a sorted copy, but that doubles
//   the memory for 'bvv' unecessarily.

void ReadsToPathsCoreY(BaseVecVec          & bvv, 
                       const size_t          K, 
                       vecKmerPath         & paths, 
                       const String        & PARCEL_DIR = "",
                       const size_t          NUM_THREADS = 1,
                       const String          CHECKPOINT = "",
                       const bool            VERBOSE = false);





// These versions of ReadsToPathsCoreY generate paths_rc and pathsdb as well.

void ReadsToPathsCoreY(BaseVecVec          & bvv, 
                       const size_t          K, 
                       vecKmerPath         & paths, 
                       vecKmerPath         & paths_rc,
                       vec<tagged_rpint>   & pathsdb,
                       const String        & PARCEL_DIR = "",
                       const size_t          NUM_THREADS = 1,
                       const String          CHECKPOINT = "",
                       const bool            VERBOSE = false);



void ReadsToPathsCoreY(BaseVecVec            & bvv, 
                       const size_t            K, 
                       vecKmerPath           & paths, 
                       vecKmerPath           & paths_rc,
                       vec<big_tagged_rpint> & pathsdb, 
                       const String          & PARCEL_DIR = "",
                       const size_t            NUM_THREADS = 1,
                       const String            CHECKPOINT = "",
                       const bool              VERBOSE = false);







#endif
