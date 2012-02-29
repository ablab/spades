//////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATCHER_COTTAGE_CORE_H
#define PATCHER_COTTAGE_CORE_H

#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"

class PatcherCottage_LogParams {

     public:

     double verbosity;
     Bool log_correct;
     Bool log_align_reads;
     Bool log_aligns;
     Bool log_aligns_all;

     PatcherCottage_LogParams( ) : verbosity(0), log_correct(False), 
          log_align_reads(False), log_aligns(False), log_aligns_all(False) { }

     PatcherCottage_LogParams( const double verbosity, const Bool log_correct,
          const Bool log_align_reads, const Bool log_aligns,
          const Bool log_aligns_all ) : verbosity(verbosity), 
          log_correct(log_correct), log_align_reads(log_align_reads),
          log_aligns(log_aligns), log_aligns_all(log_aligns_all) { }

};

void PatcherCottageCore( basevector L, basevector R, const int sep,
     const int dev, vecbasevector& reads, vecqualvector& quals, 
     vec< pair<int,int> >& pairs, String& report, 
     const PatcherCottage_LogParams& log_params, const String& data_dir, 
     vec<int>& start_stop, const int K, const int MAX_READS, const int MAX_JOINS, 
     const int MIN_OVERLAP_END, const String& ROOT, size_t p,
     const int u1, const int u2, const Bool LR_special );

#endif
