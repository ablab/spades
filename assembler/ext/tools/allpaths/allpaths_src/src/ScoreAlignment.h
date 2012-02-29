///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SCOREALIGNMENTONE
#define SCOREALIGNMENTONE

#include "Alignment.h"
#include "math/Arith.h"
#include "Basevector.h"
#include "Qualvector.h"

Float ScoreAlignment( const align& a, const basevector& rd1, 
     const qualvector& scores1, const basevector& rd2, 
     const qualvector& scores2 = qualvector(0), int start1 = 0,
     int stop1 = -1, int start2 = 0, int stop2 = -1, Bool ignore_gaps = False );

Float ScoreAlignment( Bool rd2_is_rc, const align& a, const basevector& rd1, 
     const qualvector& scores1, const basevector& rd2, 
     const qualvector& scores2 = qualvector(0), int start1 = 0,
     int stop1 = -1, int start2 = 0, int stop2 = -1, Bool ignore_gaps = False );

int ScoreAlignmentPoly( const align& a, const basevector& rd1, 
     const qualvector& scores1, const basevector& rd2, 
     const qualvector& scores2 = qualvector(0), int start1 = 0,
     int stop1 = -1, int start2 = 0, int stop2 = -1 );

int ScoreAlignmentPoly( Bool rd2_is_rc, const align& a, const basevector& rd1, 
     const qualvector& scores1, const basevector& rd2, 
     const qualvector& scores2 = qualvector(0), int start1 = 0,
     int stop1 = -1, int start2 = 0, int stop2 = -1 );

void Regap( align& a, 
	    const basevector& rd1, const qualvector& scores1,
	    const basevector& rd2, const qualvector& scores2 );

void Regap( Bool rd2_is_rc, align& a, 
	    const basevector& rd1, const qualvector& scores1, 
	    const basevector& rd2, const qualvector& scores2 );

inline 
void Regap( alignment& a, 
	    const basevector& rd1, const qualvector& scores1, 
	    const basevector& rd2, const qualvector& scores2 )
{    
  align al = align(a);
  Regap( al, rd1, scores1, rd2, scores2 );
  a.Set( packalign(al), a.Errors( ) );    
}

#endif
