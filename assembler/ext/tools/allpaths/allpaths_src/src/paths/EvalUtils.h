/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef EVAL_UTILS_H
#define EVAL_UTILS_H

#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "math/HoInterval.h"
#include "lookup/LookAlign.h"
#include "paths/AlignHyperKmerPath.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"



void SummarizeReferenceCoverage
( longlong& total_bases, longlong& total_covered, ostream& out,
  const vecbasevector& genome, const vecbasevector& genome_diploid,
  const vecbitvector& genome_amb, const vec<look_align>& aligns,
  const Bool brief = False );

void PrintGraphStatistics
( ostream& out, const HyperKmerPath& h,
  const vecbasevector& genome, const vecbasevector& genome_diploid,
  const vec<look_align>& aligns, const vec< vec<int> >& aligns_index );


// EvaluateAssembly.  Print out an evaluation of the assembly, using truth data if
// available.  Note that this will reorder the HyperKmerPath h.

void EvaluateAssembly( HyperKmerPath& h, const KmerBaseBroker* kbb,
       const String& data_dir, const String& wrun_dir, const String& sub_dir,
       const vecbasevector& genome, const vecbitvector& genome_amb, 
       const Bool DIPLOID, const Bool USE_TRUTH,
       const Bool FILTER_ALIGNS, const Bool WRITE_ALIGNS, const Bool REORDER = True,
       nbases_t MIN_TRUSTED_PATH_LEN = 0, const filename_t& HYPER = "hyper", 
       String report_suffix = "" );

void FilterAligns( const HyperKmerPath& h, vec<look_align>& aligns,
     vec< vec<int> >& aligns_index, vec<TrustedPath>& trusted_paths,
     const int MIN_LEN = 3000 );

/**
   FuncDecl: ComputeComponentSizes
   
   Compute the component sizes.  The size of each
   component here is defined to be the sum of the component's
   edges, which is not quite right but usually close enough.  To
   compensate for duplicate edges, we only take the longest edge
   between any two vertices.

   Output params:

     component_sizes - the size of each component; sorted, so you can call N50() on it.
   
 */
void ComputeComponentSizes( const HyperKmerPath& h, vec< nkmers_t >& component_sizes );

#endif
