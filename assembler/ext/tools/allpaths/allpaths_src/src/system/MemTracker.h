///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This code allows some tracking of news and deletes.

#ifndef MEMORY_TRACKER_H
#define MEMORY_TRACKER_H

// This function is registered by RunTime() as an atexit function.  It
// prints out all the memory blocks that were newed but not deleted at
// the end of the program's execution.  If memory is not being
// tracked, it has no effect.
void
check_memory_tracker_records();

// Check the bookends of all of the memory allocations we're tracking.
// If any have been corrupted, print them out, and then a list of all
// currently tracked memory allocations, which may help figure out
// where the clobbering is happening.
void
check_malloc_bookends();

// Return tracked memory currently in use
long
get_memory_usage();

// If you want to track news and deletes in some module, #define
// TRACK_MEMORY and then #include "MemTracker.h".  That will attempt
// to shunt calls to new and delete to these home-brewed versions
// below.  This allows us to check deletes against a list of newed
// objects to ensure we're not double-deleting pointers.

#ifdef TRACK_MEMORY

// We have to include these here because STL uses placement new all
// over the place which, since it takes an argument, is messed up by
// the macro that follows.

// This macro allows calls to operator new to contain debugging info.

#include <new>
#include <memory>
#include <string>

void * 
operator new ( size_t size, const char *file, const int line ) throw();

void * 
operator new[] ( size_t size, const char *file, const int line ) throw();

// This macro allows calls to operator new to contain debugging info.
#define new       new(__FILE__,__LINE__)

#endif

#endif
