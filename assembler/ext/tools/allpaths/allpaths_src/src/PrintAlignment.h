///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PRINTALIGNMENT
#define PRINTALIGNMENT

#include <fstream>

#include "Alignment.h"
#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"

void PrintBlanks( ostream& out, int n );

template<class BASEVEC>
void PrintBases( ostream& out, const BASEVEC& rd, int from, int to );

template<class BASEVEC1, class BASEVEC2>
void PrintVisualAlignment( Bool abbreviate, ostream& out, const BASEVEC1& rd1, 
     const BASEVEC2& rd2, const align& a, 
     const qualvector& scores1 = qualvector(0), 
     const qualvector& scores2 = qualvector(0), 
     int begin = 0, Bool one_frame = false, int min_score_to_abbrev = 0,
     Bool abbeviate_poor = False, float min_fract_poor = 2.0,
     Bool abbreviate_good = False, float max_fract_good = 0.05,
     Bool print_heads_and_tails = True, const Bool CtoT_special = False );

template<class BASEVEC1, class BASEVEC2>
void PrintVisualAlignment( Bool rd2_is_rc, Bool abbreviate, ostream& out, 
     const BASEVEC1& rd1, BASEVEC2 rd2, const align& a, 
     const qualvector& scores1 = qualvector(0), qualvector scores2 = qualvector(0),
     int begin = 0, Bool one_frame = false, int min_score_to_abbrev = 0,
     Bool abbeviate_poor = False, float min_fract_poor = 2.0,
     Bool abbreviate_good = False, float max_fract_good = 0.05 );

template<class BASEVEC1, class BASEVEC2>
inline void PrintVisualAlignment( Bool abbreviate, ostream& out, 
     const BASEVEC1& rd1, const BASEVEC2& rd2, const alignment& a, 
     const qualvector& scores1 = qualvector(0), 
     const qualvector& scores2 = qualvector(0),
     int begin = 0, Bool one_frame = false, int min_score_to_abbrev = 0,
     Bool abbreviate_poor = False, float min_fract_poor = 2.0 )
{    PrintVisualAlignment( abbreviate, out, rd1, rd2, align(packalign(a)),
          scores1, scores2, begin, one_frame, min_score_to_abbrev,
          abbreviate_poor, min_fract_poor );    }

template<class BASEVEC1, class BASEVEC2>
inline void PrintVisualAlignment( Bool rd2_is_rc, Bool abbreviate, ostream& out, 
     const BASEVEC1& rd1, BASEVEC2 rd2, const alignment& a, 
     const qualvector& scores1 = qualvector(0),
     qualvector scores2 = qualvector(0),
     int begin = 0, Bool one_frame = false, int min_score_to_abbrev = 0,
     Bool abbreviate_poor = False, float min_fract_poor = 2.0 )
{    PrintVisualAlignment( rd2_is_rc, abbreviate, out, rd1, rd2, 
          align(packalign(a)), scores1, scores2, begin, one_frame,
          min_score_to_abbrev, abbreviate_poor, min_fract_poor );    }

void PrintAlignment( ostream& out, const basevector& rd1, 
     const basevector& rd2, const alignment& a );
void PrintAlignment( Bool rd2_is_rc, ostream& out, const basevector& rd1, 
     basevector rd2, const alignment& a );
void PrintReadWithScores( basevector& B, qualvector& Q, ostream& out,
     const int pw = 80 );
void PrintErrorsInAlignment(ostream& out, const basevector& rd1, 
     const basevector& rd2, const alignment& a, const qualvector& scores1,
     const qualvector& scores2 );

#endif
