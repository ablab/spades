///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file LockedData.cc
 * \author tsharpe
 * \date Jul 17, 2009
 *
 * \brief
 */
// MakeDepend: library PTHREAD
#include "system/LockedData.h"
#include "system/ErrNo.h"
#include "system/System.h"

void LockedData::die( int errNo, char const* msg )
{
    ErrNo err(errNo);
    FatalErr(msg << err);
}
