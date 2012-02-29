///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef CORRECT_LONG_READS_TOOLS_H
#define CORRECT_LONG_READS_TOOLS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "PackAlign.h"
#include "math/HoInterval.h"
#include "paths/LongReadTools.h"
#include "paths/Uniseq.h"

class heuristics {

     public:

     int L;
     int min_overlap_to_see_other;
     int max_offset_diff_for_other;
     int min_rdist_for_other;
     int delta_r_max;
     double max_delta_ratio;
     int delta_sub;
     int min_spread;
     double min_ratio_to_kill;
     int max_paths;
     int max_iterations;
     int min_safe_cn1;
     double min_win_ratio;
     double min_votes;
     double min_votes_to_protect;

};

void FollowVotes( int i, int low, int high,
     const int np, const vec< vec<int> >& ppaths, const vec<double>& votes,
     const vec< pair<int,int> >& verts, const vecbasevector& unibases,
     const int K, const heuristics& heur, ostringstream& hout,
     vec< vec<int> >& answer, vec<double>& answer_votes );

void ComputeMatches( const vec<int>& p, const vecbasevector& unibases,
     const vec< pair<int,int> >& verts, const vec<int>& unis,
     const vec< vec< pair<int,int> > >& umatches, ostream& out,
     const Bool PRINT_MATCHES, vec< pair<int,int> >& matches );

void AlignReadToPath( const basevector& r, const int id, const Bool fw, 
     const vec<int>& p, const int path_id, const vecbasevector& unibases,
     const vec< pair<int,int> >& verts, const vec<int>& unis,
     const vec< vec< pair<int,int> > >& umatches, const vec< vec<int> >& ppaths,
     vec<Bool>& ppath_seen, const int L, const Bool PRINT_MATCHES,
     const Bool DUMP_LOCAL, const Bool PRINT_LM1, double& ERRS, int& LMATCHES,
     triple<int,int,String>& EREPORTS, int& LASTPOS, int& LAST, Bool& SKIPPED,
     map< pair<ho_interval,basevector>, pair<align,double> >& A );

void Phase2( const int u, const vecbasevector& unibases, 
     const vec< vec< pair<int,int> > >& Ulocs, const vec<int>& to_rc,
     const vecbasevector& longreads, const int K, const vec< vec<int> >& hits, 
     const vec<uniseq>& UNISEQ, const vec< vec<int> >& UNISEQ_ID, ostream& hout, 
     const Bool DOT2, 
     const Bool PRINT_MATCHES, const Bool PRINT_LM1, const Bool DUMP_LOCAL, 
     const Bool PRINT_DISCARDS, const Bool VALIDATE1, const String& data_dir, 
     const vecbasevector& genome2, const vec< vec< pair<int,int> > >& Glocs, 
     const int LG, const heuristics& heur, vecbasevector& new_stuff,
     vec< vec<int> >& right_exts,
     const Bool NEW_FILTER, const Bool ANNOUNCE, const int SUPER_VERBOSITY,
     const Bool PRINT_BASES2, const int MIN_SPREAD2 );

void Validate( vecbasevector all, const vecbasevector& unibases,
     const vecbasevector& genome, const int min_mult, 
     const Bool PRINT_MISSING_KMERS );

template<int K> void BuildPatches( const vecbasevector& unibases, 
     const vec< vec<int> >& nexts, const vec<Bool>& cn_plus, const int L, 
     const vec< vec< pair<int,int> > >& Ulocs, vec<GapPatcher>& patchers, 
     const Bool CORRECT_PATCHES, const Bool CORRECT_PATCHES_NEW, 
     const Bool CORRECT_PATCH_VERBOSE, const Bool ANNOUNCE_PATCH_CORRECTION,
     const vecbasevector& genome2, const int PATCH_VERBOSITY, 
     const Bool VALIDATE_PATCHES, const int PATCH_CORRECT_VERBOSITY, const int LG, 
     const vec< vec< pair<int,int> > >& Glocs, const String& data_dir, 
     const String& run_dir, vec<basevector>& bpatches, const uint n_patchers_min,
     const vecbasevector& fbases, const vecqualvector& fquals,
     const PairsManager& fpairs, const vec< kmer<20> >& fheads, 
     const vec<int64_t>& fids );

#endif
