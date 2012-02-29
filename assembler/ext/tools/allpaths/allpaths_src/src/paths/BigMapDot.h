///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BIG_MAP_DOT_H
#define BIG_MAP_DOT_H

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/BigMapTools.h"
#include "paths/Uniseq.h"

void BigMapDot( const String& DOT, const snark& S, const int K2, 
     const Bool VALIDATE, const vecbasevector& genome, 
     const vec< vec<placementy> >& Ulocs2, const vec<String>& legends_to_show,
     const Bool circo );

#endif
