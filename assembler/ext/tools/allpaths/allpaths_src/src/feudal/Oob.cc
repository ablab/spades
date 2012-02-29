///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Oob.cc
 * \author tsharpe
 * \date Jul 2, 2009
 *
 * \brief Just a function to call FatalErr on out-of-bounds access.
 */
#include "feudal/Oob.h"
#include "system/System.h"
#include <cstdlib>

void OutOfBoundsReporter::oob( char const* className, size_t idx, size_t siz )
{
    FatalErr(className << " access out of bounds.  Tried to access " << idx << " when there were only " << siz << " elements.");
}
