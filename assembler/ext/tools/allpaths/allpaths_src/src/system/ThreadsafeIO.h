///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file ThreadsafeIO.h
 * \author tsharpe
 * \date May 14, 2009
 *
 * \brief Provides threadsafe wrappers to control access to an ostream by multiple threads.
 *
 * To use the static instances that wrap cout and cerr, replace cout in your code with ThreadsafeOStreamFetcher::cout()
 * and cerr with ThreadsafeOStreamFetcher::cerr().  You can define a couple of macros to make this easier,
 * if you like, e.g., #define cout ThreadsafeOStreamFetcher::cout()
 */
#ifndef SYSTEM_THREADSAFEIO_H_
#define SYSTEM_THREADSAFEIO_H_

#include <ostream>
#include <map>

/// Streambuf which locks a mutex and then writes to a wrapped ostream on overflow and sync.
class ThreadsafeStreambuf : public std::streambuf
{
public:
    ThreadsafeStreambuf( std::ostream& os, pthread_mutex_t* pLock )
    : mOS(os), mpLock(pLock)
    { setp(mBuf,mBuf+sizeof(mBuf)-1); }

    ~ThreadsafeStreambuf()
    { sync(); }

private:
    ThreadsafeStreambuf( ThreadsafeStreambuf const& ); // undefined -- no copying
    ThreadsafeStreambuf& operator=( ThreadsafeStreambuf const& ); // undefined -- no copying

    int_type overflow( int_type ch );
    int sync();

    std::ostream& mOS;
    pthread_mutex_t* mpLock;
    char mBuf[4096];
};

/// An ostream that wraps another.
/// If each thread uses a separate one of these, then access to the wrapped ostream becomes threadsafe.
class ThreadsafeOStream : public std::ostream
{
public:
    ThreadsafeOStream( std::ostream& os, pthread_mutex_t* pLock )
    : std::ostream(&mSB), mSB(os,pLock)
    {}

private:
    ThreadsafeOStream( ThreadsafeOStream const& ); // undefined -- no copying
    ThreadsafeOStream& operator=( ThreadsafeOStream const& ); // undefined -- no copying

    ThreadsafeStreambuf mSB;
};

/// A class which looks up (creating, if necessary) a ThreadsafeOStream for the current thread.
/// Each instance of this class is specific to a particular ostream being accessed by multiple threads.
/// <br>
/// We make a static instance to wrap cout and cerr, for convenience.  You could make others to control,
/// for example, access to a file that is being written by multiple threads.
/// <br>
/// Special note on cleanup, and the destroy function:<br>
/// Officially, you have to call destroy on any ThreadsafeOStreamFetcher you've used in a thread, ideally, somewhere near
/// the end of your thread function.  (This is necessary because there's no good, reliable hook to capture thread
/// exit.  We could rewrite to use boost someday to automate this.  Sorry.)<br>
/// See the example in the file WorklistTest, where we make a nice macro so that it looks like we're just using cout
/// normally, and the Processor's destructor does the correct, official cleanup.<br>
/// Unofficially, it probably doesn't matter, unless you're terminating lots of threads during the course
/// of your process's execution.  (Which is unusual:  usually you make a single thread pool, and all the threads live
/// pretty much right up to process exit.)  In this case you'll temporarily leak the ThreadsafeOStreams you've created
/// for threads that no longer exist.  The thread-specific streams will eventually get cleaned up when the Fetcher dies,
/// but that might not be prompt, especially in the case of the gOUT and gERR Fetchers, which live until program termination.<br>
/// Note that this isn't typically a problem for a ThreadsafeOStreamFetcher you make to control access to your own output file:
/// you'd typically create the ThreadsafeOStreamFetcher, create a pool of threads and let them do their thing, and then
/// destroy the Fetcher when the work is all done, promptly reclaiming resources.  It's just the cout and cerr wrappers
/// that have an unusually long lifetime that might cause problems.
class ThreadsafeOStreamFetcher
{
public:
    static std::ostream& cout() { return gCOUT.get(); }
    static void coutDone() { gCOUT.destroy(); }
    static std::ostream& cerr() { return gCERR.get(); }
    static void cerrDone() { gCERR.destroy(); }

    ThreadsafeOStreamFetcher( std::ostream& os );
    ~ThreadsafeOStreamFetcher();

    std::ostream& get(); // get a thread-specific ostream
    void destroy(); // destroy the thread-specific ostream

private:
    std::ostream& mOS;
    pthread_mutex_t mStreamLock;
    pthread_mutex_t mMapLock;
    std::map<pthread_t,ThreadsafeOStream*> mMap;

    static ThreadsafeOStreamFetcher gCOUT;
    static ThreadsafeOStreamFetcher gCERR;
};

#endif /* SYSTEM_THREADSAFEIO_H_ */
