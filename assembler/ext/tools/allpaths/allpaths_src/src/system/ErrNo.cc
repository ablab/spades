///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file ErrNo.cc
 * \author tsharpe
 * \date Dec 13, 2011
 *
 * \brief
 */
#include "system/ErrNo.h"
#include <cerrno>
#include <cstring>
#include <sstream>

ErrNo::ErrNo()
: mErrNo(errno)
{
}

// we don't know whether we're going to get the GNU version or the XOPEN version
// of strerror_r, so we're going to code so that either will work
std::string ErrNo::text() const
{
    char buf[8192];
    long val = (long)strerror_r(mErrNo,buf,sizeof(buf));

    char const* msg = "Unknown error"; // XOPEN version failed
    if ( val == 0 ) msg = buf; // XOPEN version worked
    else if ( val != -1 ) msg = (char*)val; // GNU version worked

    std::ostringstream os;
    os << ": " << msg << " [errno=" << mErrNo << "].";
    return os.str();
}
