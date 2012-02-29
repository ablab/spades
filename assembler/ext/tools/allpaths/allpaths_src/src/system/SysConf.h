///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file SysConf.h
 * \author tsharpe
 * \date Apr 15, 2010
 *
 * \brief
 */
#ifndef SYSTEM_SYSCONF_H_
#define SYSTEM_SYSCONF_H_

#include <cstddef>

/// The number of bytes per page. (Sometimes things are measured in pages.)
size_t pageSize();

/// Number of bytes of memory on this machine.
size_t physicalMemory();

/// Number of bytes of memory on this machine that are not currently in use
/// by user processes.
size_t availablePhysicalMemory();

/// The number of CPUs available.
size_t processorsOnline();

/// The maximum number of threads to use, given all constraints.
int configNumThreads( int numThreads );

/// The number of ticks per second.  (Sometimes things are measured in ticks.)
size_t clockTicksPerSecond();

/// The maximum length of a host name (for gethostname).
size_t maxHostNameLen();

#endif /* SYSTEM_SYSCONF_H_ */
