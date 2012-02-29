/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file LookupTab.cc
 * \author tsharpe
 * \date Jun 15, 2009
 *
 * \brief Instantiates Location vectors.
 *
 *
 */
#include "lookup/LookupTab.h"
#include "feudal/OuterVecDefs.h"
#include "feudal/SmallVecDefs.h"

template class SmallVec< Location, MempoolAllocator<Location> >;
template class OuterVec<LocationVec>;
