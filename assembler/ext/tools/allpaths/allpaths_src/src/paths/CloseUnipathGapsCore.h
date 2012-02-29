/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef ALIGN_TO_UNIGAPS_CORE
#define ALIGN_TO_UNIGAPS_CORE

#include "Basevector.h"
#include "Intvector.h"
#include "Qualvector.h"
#include "STLExtensions.h"
#include "String.h"
#include "Vec.h"

// Fill the seeds vector.
void FindUnipathGapSeeds(const vecbasevector & bases, 
                         const vecbasevector & unibases,
                         const size_t K, 
                         const size_t n_threads, 
                         const String & work_dir, 
                         VecULongVec * seeds);


void CloseUnipathGapsCore( const vecbasevector & bases, 
                            const vecqualvector & quals,
                            const vecbasevector & unibases, 
                            const vec< vec<int> > & nexts, 
                            const vec<int> & to_rc, 
                            const size_t K, 
                            const size_t UNIBASES_K,
                            const VecULongVec & seeds,
                            const VecULongVec & seeds_rc,
                            vec< triple<int,int,longlong> > & extenders, 
                            const int VERBOSITY, 
                            ostream & log);

#endif
