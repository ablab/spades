///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Exit.cc
 * \author tsharpe
 * \date Feb 6, 2009
 *
 * \brief A replacement for the ::exit(int) function.
 */
// MakeDepend: library PTHREAD
#include "system/Exit.h"
#include <stdlib.h>
#include <pthread.h>

namespace
{

class MainThread
{
public:
    MainThread() : mThread(pthread_self()) {}
    pthread_t val() { return mThread; }
private:
    pthread_t mThread;
};

// since this is created at static initialization time, its val
// will be that of the main thread.
MainThread gMainThread;

} // end of anonymous namespace


void CRD::exit( int status )
{
    if ( !pthread_equal(gMainThread.val(),pthread_self()) )
    {
        /// TODO: It would be awfully nice to produce a traceback of the thread.
        pthread_exit( reinterpret_cast<void*>(status) );
    }
    if ( !status )
        ::exit(0);
    abort();
}
