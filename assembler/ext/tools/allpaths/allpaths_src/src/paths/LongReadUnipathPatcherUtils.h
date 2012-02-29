//////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__LONG_READ_UNIPATH_PATCHER_UTILS__H
#define PATHS__LONG_READ_UNIPATH_PATCHER_UTILS__H

#include "Basevector.h"
#include "paths/LongReadTools.h"

void BuildGLocs( const int LG,
		 const vecbasevector &genome,
		 vec< vec< pair<int,int> > > &Glocs );

void AlignToGenome( const int LG,
		    const bool COMPUTE_TRUE_READ_LENGTH,
		    const basevector &r,
		    const vecbasevector &genome,
		    const vec< vec< pair<int,int> > > &Glocs,
		    vec<align_data> &adata,
		    ostream &out );

#endif
