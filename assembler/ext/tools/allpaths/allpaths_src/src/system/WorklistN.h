/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2010) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file WorklistN.h
 * \author tsharpe
 * \date Apr. 14, 2010
 *
 * \brief Use multiple threads of execution to process a number of workitems.
 * This is just like the more general Worklist, but the workitems are simply
 * numbered, rather than existing as Workitem objects in a list.
 */
#ifndef SYSTEM_WORKLISTN_H_
#define SYSTEM_WORKLISTN_H_

#include "system/Assert.h"
#include "system/SysConf.h"
#include "system/WorklistUtils.h"
#include <algorithm>
#include <cstddef>
#include <list>
#include <time.h>

/// See the documentation for Worklist.  This works the same way.
template<class Processor>
class WorklistN
{
public:
    struct ProgressMonitor
    {
        ProgressMonitor() : mWorkitem(~0ul), mTime(0) {}
        size_t mWorkitem;
        time_t mTime;
    };

    /// Create a Worklist that will be processed by nThreads of execution.
    /// If nThreads is 0, no threads will be created, and the entire
    /// Worklist will be processed by the Worklist destructor.  The default
    /// value for nThreads is the number of processors less one (we assume that
    /// the main thread will continue to do work).
    /// Passing a stack size of 0 (it's default value) will cause threads to be
    /// created with a stack as large as main thread's stack (unless the main
    /// thread's stack size is unlimited, in which case a default value of 8Mb
    /// will be used for each thread's stack size).
    WorklistN( Processor const& p,
               size_t nWorkitems,
               size_t nThreads = processorsOnline()-1,
               size_t threadStackSize = 0 )
    : mCurWorkitem(0ul), mNWorkitems(nWorkitems),
      mNThreads(calcThreads(nThreads,nWorkitems)), mProcessor(p),
      mProgressMonitors(new ProgressMonitor[mNThreads+1]),
      mTP(mNThreads,threadStackSize,threadFunc,this)
    { QueueStateManipulator qsm(mQS); qsm.setSize(nWorkitems); }

    /// This will wait until all work is done, and all threads have
    /// been terminated.  That might be a long time!
    ~WorklistN()
    { waitForDone();
      delete [] mProgressMonitors; }

    /// Start up a separate thread to monitor thread death.
    /// This thread wakes up every howOftenToCheckInSecs seconds, and replaces
    /// any threads that have died with a new thread.
    /// If you call this more than once, subsequent calls merely adjust the
    /// check interval.
    /// If you call this with howOftenToCheckInSecs == 0, the monitor exits.
    /// If reportRestarts is true, an error message is written to cerr to inform
    /// you of a thread restart.
    void threadDeathMonitor( long howOftenToCheckInSecs,
                             bool reportRestarts = true )
    { mTP.threadDeathMonitor(howOftenToCheckInSecs,reportRestarts); }

    /// Return the number of unprocessed workitems.
    size_t size()
    { QueueStateManipulator qsm(mQS);
      return qsm.getSize(); }

    /// Note:  you don't have to call this -- the destructor will automatically
    /// wait until all the work is complete.  You may call this for
    /// entertainment value, if you'd like.
    void waitForDone()
    { mTP.threadDeathMonitor(0); // shut down the monitor.
      QueueStateManipulator(mQS).setDone(); // tell threads to quit
      doWork(true);
      mTP.shutdown(); }

    /// Clear the worklist of all unprocessed items.
    /// If you need to quit early, this is the only way to do it.
    /// The destructor is still going to wait for in-progress work to complete.
    void clear()
    { QueueStateManipulator qsm(mQS);
      mCurWorkitem = mNWorkitems;
      qsm.setSize(0); }

    /// How many threads are there to process the work.
    size_t getNThreads()
    { return mTP.getNThreads(); }

    /// What happening with a thread.
    /// ThreadIdx runs from 0 to getNThreads()-1.
    ProgressMonitor getProgress( size_t threadIdx )
    { QueueStateManipulator qsm(mQS);
      return mProgressMonitors[threadIdx]; }

private:
    // unimplemented -- no copying
    WorklistN( WorklistN const& );
    WorklistN& operator=( WorklistN const& );

    void doWork( bool isMainThread = false );
    static void* threadFunc( void* );

    static size_t calcThreads( size_t nThreads, size_t nWorkitems )
    { return nWorkitems ? std::min(nThreads+1,nWorkitems)-1 : 0; }

    size_t mCurWorkitem;
    size_t mNWorkitems;
    size_t mNThreads;
    QueueState mQS;
    Processor mProcessor;
    ProgressMonitor* mProgressMonitors;
    ThreadPool mTP; // this needs to be the last member because we don't
                    // want to start threads until everything else is ready
};

template<class Processor>
void WorklistN<Processor>::doWork( bool isMainThread )
{
    size_t idx = mTP.getNThreads();
    if ( !isMainThread )
        idx = mTP.findThreadIndex(pthread_self());
    ProgressMonitor& progressMonitor = mProgressMonitors[idx];
    Processor processor(mProcessor);

    while ( true )
    {
        if ( true ) // empty block just to control lifetime of qsm
        {
            QueueStateManipulator qsm(mQS);
            qsm.waitForWork();
            if ( !qsm.getSize() ) // wait for work only returns when there's work, or when the done flag is set
                break;            // so break if there's no work (because that means we're all done)

            time(&progressMonitor.mTime);
            progressMonitor.mWorkitem = mCurWorkitem++;
            qsm.decSize();
        }
        processor(progressMonitor.mWorkitem);
    }
}

template<class Processor>
void* WorklistN<Processor>::threadFunc( void* ptr )
{
    reinterpret_cast<WorklistN*>(ptr)->doWork();
    return 0;
}

#endif /* SYSTEM_WORKLISTN_H_ */
