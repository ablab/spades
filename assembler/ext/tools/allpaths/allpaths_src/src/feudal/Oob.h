///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Oob.h
 * \author tsharpe
 * \date Jul 2, 2009
 *
 * \brief Just a function to call FatalErr on out-of-bounds access.
 * The only purpose is to avoid having to include the kitchen-sink via system/System.h
 * in classes that would be better off being spare.
 */

#ifndef FEUDAL_OOB_H_
#define FEUDAL_OOB_H_

#include <cstddef>

struct OutOfBoundsReporter
{
    static void oob( char const* className, size_t idx, size_t siz );
};

#endif /* FEUDAL_OOB_H_ */
