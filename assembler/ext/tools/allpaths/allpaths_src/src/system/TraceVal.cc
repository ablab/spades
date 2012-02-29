///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <iostream>
#include "system/TraceVal.h"

TraceValCommon::timestamp_t TraceValCommon::nextTimeStamp_ = 1;
Bool TraceValCommon::tracingCopiesStopped_ = False;
TraceValCommon::timestamp_t TraceValCommon::stopAtStamp_ = 0;

template class TraceVal<int>;
