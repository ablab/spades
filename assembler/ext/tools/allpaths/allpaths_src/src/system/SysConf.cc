///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file SysConf.cc
 * \author tsharpe
 * \date Apr 15, 2010
 *
 * \brief
 */
#include "system/SysConf.h"
#include "system/Exit.h"
#include <climits>
#include <cstdlib>
#include <iostream>

namespace
{

void sysconfErr( char const* attr )
{
    std::cout << "Sysconf unable to determine value of " << attr << std::endl;
    CRD::exit(1);
}

size_t gPagSiz;
size_t gPhysPages;
size_t gProcsOnline;
size_t gClockTicksPerSec;
size_t gMaxHostNameLen;

}

size_t pageSize()
{
    if ( !gPagSiz )
    {
        long result = sysconf(_SC_PAGESIZE);
        if ( result == -1 ) sysconfErr("_SC_PAGESIZE");
        gPagSiz = result;
    }
    return gPagSiz;
}

size_t physicalMemory()
{
    if ( !gPhysPages )
    {
        long result = sysconf(_SC_PHYS_PAGES);
        if ( result == -1 ) sysconfErr("_SC_PHYS_PAGES");
        gPhysPages = result;
    }
    return gPhysPages*pageSize();
}

size_t availablePhysicalMemory()
{
  long result = sysconf(_SC_AVPHYS_PAGES);
  if ( result == -1 ) sysconfErr("_SC_AVPHYS_PAGES");
  return result * pageSize();
}

size_t processorsOnline()
{
    if ( !gProcsOnline )
    {
        long result = sysconf(_SC_NPROCESSORS_ONLN);
        if ( result == -1 ) sysconfErr("_SC_NPROCESSORS_ONLN");
        if ( result < 1 || result > INT_MAX )
        {
            std::cout << "Sysconf's value for the number of processors online ("
                      << result << ") doesn't make sense." << std::endl;
            CRD::exit(1);
        }
        gProcsOnline = result;
    }
    return gProcsOnline;
}

// If the argument is silly (meaning not positive), then we check the
// environment variable OMP_THREAD_LIMIT for inspiration.  if that's not set
// we go with the number of processors online.  we do a final check to make
// sure that, whereever we got the number, that it isn't set to more than the
// number of processors online.
int configNumThreads( int numThreads )
{
    int procsOnline = processorsOnline();
    if ( numThreads < 1 )
    {
        char const* ompEnv = getenv("OMP_THREAD_LIMIT");
        if ( !ompEnv )
            numThreads = procsOnline;
        else
        {
            char* end;
            long ompVal = strtol(ompEnv,&end,10);
            if ( *end || ompVal < 1 || ompVal > INT_MAX )
            {
                std::cout << "Environment variable OMP_THREAD_LIMIT ("
                          << ompEnv
                          << " cannot be parsed as a positive integer."
                          << std::endl;
                CRD::exit(1);
            }
            numThreads = ompVal;
        }
    }
    if ( numThreads > procsOnline )
        numThreads = procsOnline;

    return numThreads;
}

size_t clockTicksPerSecond()
{
    if ( !gClockTicksPerSec )
    {
        long result = sysconf(_SC_CLK_TCK);
        if ( result == -1 ) sysconfErr("_SC_CLK_TCK");
        gClockTicksPerSec = result;
    }
    return gClockTicksPerSec;
}

size_t maxHostNameLen()
{
    if ( !gMaxHostNameLen )
    {
        long result = sysconf(_SC_HOST_NAME_MAX);
        if ( result == -1 ) sysconfErr("_SC_HOST_NAME_MAX");
        gMaxHostNameLen = result;
    }
    return gMaxHostNameLen;
}
