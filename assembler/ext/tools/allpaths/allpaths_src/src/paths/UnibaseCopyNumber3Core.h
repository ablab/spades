///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// UnibaseCopyNumber3Core: predict copy numbers.  Warning: this creates files
// in a way which it shouldn't.  This feature will presumably be eliminated in
// a subsequent version.

#ifndef UNIBASE_COPY_NUMBER3_CORE
#define UNIBASE_COPY_NUMBER3_CORE

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/PdfEntry.h"

const double THRESH_DEFAULT = 0.01;
const double ERR_RATE_DEFAULT = 0.0;
const Bool CORRECT_BIAS_DEFAULT = True;

void UnibaseCopyNumber3Core( 
     // INPUTS:
     const int K, const vecbasevector& reads, const vecbasevector& unibases, 
     const int NUM_THREADS, const String& reads_head, const String& unibases_head,
     const int PLOIDY, 
     // LOGGING CONTROL:
     const Bool VERBOSE, ostream& logout,
     // OUTPUT:
     VecPdfEntryVec& cn_pdfs,
     vec<double>& cn_raw,
     vec<int64_t>& n_kmer_hits_raw,
     // HEURISTIC PARAMETERS:
     const double THRESH = THRESH_DEFAULT, const double ERR_RATE = ERR_RATE_DEFAULT,
     const Bool CORRECT_BIAS = CORRECT_BIAS_DEFAULT
     );

void ComputeProbs(
     // INPUTS:
     const int K,
     const int PLOIDY,
     double occCnPloidy,
     double THRESH,
     double ERR_RATE,
     const vecbasevector& unibases,
     const vec<int64_t>& n_kmer_hits,
     // OUTPUTS:
     VecPdfEntryVec& cn_pdfs,
     vec<double>& cn_raw,
     vec<longlong>& cnHisto,
     // LOGGING:
     const Bool VERBOSE,
     ostream& logout );

#endif
