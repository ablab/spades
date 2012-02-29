///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file ThreadsafeIO.cc
 * \author tsharpe
 * \date May 14, 2009
 *
 * \brief Redefine cout, cerr in a threadsafe manner.
 */
#include <iostream>
#include "system/ThreadsafeIO.h"

ThreadsafeOStreamFetcher ThreadsafeOStreamFetcher::gCOUT(std::cout);
ThreadsafeOStreamFetcher ThreadsafeOStreamFetcher::gCERR(std::cerr);

typedef std::map<pthread_t,ThreadsafeOStream*>::iterator MapItr;

ThreadsafeOStreamFetcher::ThreadsafeOStreamFetcher( std::ostream& os )
: mOS(os)
{
    pthread_mutex_init(&mStreamLock, 0);
    pthread_mutex_init(&mMapLock, 0);
}

ThreadsafeOStreamFetcher::~ThreadsafeOStreamFetcher()
{
    MapItr end(mMap.end());
    for ( MapItr itr(mMap.begin()); itr != end; ++itr )
        delete itr->second;
    pthread_mutex_destroy(&mMapLock);
    pthread_mutex_destroy(&mStreamLock);
}

std::ostream& ThreadsafeOStreamFetcher::get()
{
    pthread_mutex_lock(&mMapLock);

    ThreadsafeOStream* pResult = mMap[pthread_self()];
    if ( !pResult )
    {
        pResult = new ThreadsafeOStream(mOS,&mStreamLock);
        mMap[pthread_self()] = pResult;
    }

    pthread_mutex_unlock(&mMapLock);

    return *pResult;
}

void ThreadsafeOStreamFetcher::destroy()
{
    pthread_mutex_lock(&mMapLock);

    MapItr pos(mMap.find(pthread_self()));
    if ( pos != mMap.end() )
    {
        delete pos->second;
        mMap.erase(pos);
    }

    pthread_mutex_unlock(&mMapLock);
}

ThreadsafeStreambuf::int_type ThreadsafeStreambuf::overflow( int_type ch )
{
    if ( ch != traits_type::eof() )
    {
        *pptr() = ch;
        pbump(1);
    }
    return sync() ? traits_type::eof() : ch;
}

int ThreadsafeStreambuf::sync()
{
    pthread_mutex_lock(mpLock);

    char* buf = pbase();
    mOS.write(buf,pptr() - buf);
    setp(buf,epptr());

    pthread_mutex_unlock(mpLock);
    return mOS.fail(); // i.e., 0 if we wrote everything, 1 if we didn't.
}
