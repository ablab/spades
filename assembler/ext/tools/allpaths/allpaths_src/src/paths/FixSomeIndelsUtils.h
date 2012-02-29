///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__FIX_SOME_INDELS_UTILS_H
#define PATHS__FIX_SOME_INDELS_UTILS_H

#include "Basevector.h"
#include "STLExtensions.h"
#include "String.h"
#include "feudal/QualNibbleVec.h"
#include "graph/Digraph.h"
#include "paths/ReadLoc.h"
#include "polymorphism/DumbCall.h"

inline Bool cmp_pos( const triple<int64_t,int,int>& x1,
		     const triple<int64_t,int,int>& x2 )
{    if (x1.second < x2.second) return True;
     if (x1.second > x2.second) return False;
     if (x1.third < x2.third) return True;
     if (x1.third > x2.third) return False;
     return x1.first < x2.first;    }

// Select reads over a specified interval of a given contig.

void SelectReads( const Bool VERBOSE,
		  const Bool USE_QUALS,
		  const int rep_start,
		  const int rep_stop,
		  const vec< triple<int64_t,int,int> > &RALIGNS_this,
		  const vec< triple<int64_t,int,int> > &JRALIGNS_this,
		  const String &run_dir,
		  const String &TIGS,			  
		  const BaseVecVec &bases,
		  const BaseVecVec &jbases,
		  const VecQualNibbleVec &quals,
		  const VecQualNibbleVec &jquals,
		  vecbasevector &reads,
		  VecQualNibbleVec &readsq,
		  ostream *p_rout );

// FindBestHomeForReads.  Consider reads R (with quality scores Q).
// Look for gap-free alignments to "contiglet" c, such that the
// central nc2 bases of c (starting at position nc1) are covered and
// to identity at least min_core_identity.  Find the minimum number of
// errors amongst all such alignments, then score each by the sum of
// the read quality scores at mismatches.

void FindBestHomeForReads( const Bool VERBOSE,
			   const String& title,
			   const vecbasevector& reads, 
			   const VecQualNibbleVec& quals,
			   const FastaVec& c,
			   const int nc1,
			   const int nc2, 
			   const double min_core_identity,
			   vec< vec<int> >& ERRS,
			   ostream *p_rout );

// Log informations about flaky regions (unexplored Jaffe code).

void ShowFlakyRegions( const int n,
		       const int tig,
		       const double min_ref_frac,
		       const bvec &TIG,
		       const vecbvec &frag_reads,
		       const vec<dumbcall> &calls,
		       const vec<read_loc> &locs,
		       ostream &rout );


// Below are the duplicate functions created for FixAssemblyBaseErrorsa. There are some modifications
// from the original version for this situation.

void FindBestHomeForReads2( const Bool VERBOSE,
			   const String& title,
			   const vecbasevector& reads, 
			   const VecQualNibbleVec& quals,
			   const FastaVec& c,
			   const int nc1,
			   const int nc2, 
			   const double min_core_identity,
			   vec< vec<int> >& ERRS,
			   ostream *p_rout );

void SelectReads2( const Bool VERBOSE,
		  const int rep_start,
		  const int rep_stop,
		  const vec< triple<int64_t,int,int> > &RALIGNS_this,
		  const BaseVecVec &bases,
		  const VecQualNibbleVec &quals,
		  vecbasevector &reads,
		  VecQualNibbleVec &readsq,
		  ostream *p_rout );


#endif
