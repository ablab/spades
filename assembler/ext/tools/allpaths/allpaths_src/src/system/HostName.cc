///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file HostName.cc
 * \author tsharpe
 * \date Jan 9, 2012
 *
 * \brief
 */
#include "system/HostName.h"
#include "system/SysConf.h"
#include "system/System.h"

std::string getHostName()
{
    size_t len = maxHostNameLen()+1; // 1 extra for null

    char* buf = new char[len+1]; // 1 extra beyond that
    if ( gethostname(buf,len) == -1 )
        FatalErr("Unable to determine host name.");
    buf[len] = 0; // make completely, absolutely sure we're null terminated

    std::string result(buf);
    delete [] buf;
    return result;
}
