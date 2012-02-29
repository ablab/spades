///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SMITHWATFREE
#define SMITHWATFREE

#include "Alignment.h"
#include "Basevector.h"
#include "CoreTools.h"
#include "Fastavector.h"

unsigned int SmithWatFree( const basevector& S, const basevector& T, 
			   int& best_loc, alignment& a,
			   bool penalize_left_gap = false, 
			   bool penalize_right_gap = false,
                           unsigned int mismatch_penalty = 2, 
                           unsigned int gap_penalty = 3,
                           unsigned int outer_gap_penalty = 3 );

unsigned int SmithWatFree( const fastavector& S, const fastavector& T, 
			   int& best_loc, alignment& a,
			   bool penalize_left_gap = false, 
			   bool penalize_right_gap = false,
                           unsigned int mismatch_penalty = 2, 
                           unsigned int gap_penalty = 3,
                           unsigned int outer_gap_penalty = 3 );

#endif
