/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
#include "paths/simulation/Placement.h"
#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"

template class SmallVec< placement, MempoolAllocator<placement> >;
template class OuterVec<PlacementVec>;
