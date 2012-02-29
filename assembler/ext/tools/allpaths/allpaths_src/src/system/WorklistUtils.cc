///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file WorklistUtils.cc
 * \author tsharpe
 * \date Feb 4, 2009
 *
 * \brief A couple of utility classes for the Worklist.
 */
// MakeDepend: library PTHREAD
#include "system/System.h"
#include "system/WorklistUtils.h"
#include "system/Assert.h"
#include "system/ErrNo.h"
#include "system/SysConf.h"
#include <ctime>
#include <csignal>
#include <cstring>
#include <iostream>
#include <unistd.h>
#include <sys/resource.h>

WorklistParameterizer::WorklistParameterizer( size_t nUnits,
                                              size_t bytesPerUnit,
                                              size_t minBatchSize,
                                              unsigned nThreads )
{
    if ( nUnits < minBatchSize )
        minBatchSize = nUnits;

    if ( !nThreads || nThreads > processorsOnline() )
        nThreads = processorsOnline();

    size_t maxSimultaneous = .9*MemAvailable()/bytesPerUnit;
    if ( maxSimultaneous < minBatchSize )
        FatalErr("Not enough memory.");

    size_t maxThreads = maxSimultaneous/minBatchSize;
    if ( nThreads > maxThreads )
        nThreads = maxThreads;

    size_t maxBatchSize = maxSimultaneous/nThreads;
    if ( maxBatchSize > nUnits )
        maxBatchSize = nUnits;

    size_t minBatches = (nUnits+maxBatchSize-1)/maxBatchSize;
    size_t maxBatches = nUnits / minBatchSize;

    if ( maxBatches <= nThreads )
        nThreads = minBatches = maxBatches;
    else
    {
        // round up minBatches to an even divisor of nThreads, subject to
        // remaining less than maxBatches
        size_t nCycles = (minBatches+nThreads-1)/nThreads;
        size_t batches = nCycles * nThreads;
        if ( batches <= maxBatches )
            minBatches = batches;
        else
            minBatches = maxBatches;
    }
    mNBatches = minBatches;
    mNThreads = nThreads;
    mBatchSize = (nUnits+minBatches-1)/minBatches;

    ForceAssertLe(mNThreads*mBatchSize,maxSimultaneous);
}

QueueStateManipulator::~QueueStateManipulator()
{
    size_t size = mQS.getSize();
    bool done = mQS.isDone();
    mQS.unlock();

    if ( !mDone && done )
    {
        mQS.mCondNotEmpty.broadcast();
    }
    else if ( size > mSize )
    {
        if ( size-mSize == 1 )
            mQS.mCondNotEmpty.signal();
        else
            mQS.mCondNotEmpty.broadcast();
    }

    if ( mSize && !size )
    {
        mQS.mCondEmpty.broadcast();
    }
}

ThreadPool::ThreadPool( size_t nThreads, size_t threadStackSize,
                        void* (*threadFunc)(void*), void* threadFuncArg )
: mThreadFunc(threadFunc), mThreadFuncArg(threadFuncArg),
  mNThreads(validateNThreads(nThreads)), mThreads(new pthread_t[nThreads+1]),
  mCondExit(*this), mMonitorInterval(0), mReportRestarts(false)
{
    checkThreadOp(pthread_attr_init(&mAttr),
                  "ThreadAttr initialization failed: ");

    if ( !threadStackSize )
    {
        struct rlimit rl;
        if ( getrlimit(RLIMIT_STACK,&rl) )
        {
            ErrNo err;
            FatalErr("Unable to determine current stack size" << err);
        }

        if ( rl.rlim_cur == RLIM_INFINITY )
            threadStackSize = DEFAULT_THREAD_STACK_SIZE;
        else
            threadStackSize = rl.rlim_cur;
    }
    checkThreadOp(pthread_attr_setstacksize(&mAttr,threadStackSize),
                  "Unable to set stack size for threads: ");

    pthread_t* pThread = mThreads + mNThreads;
    while ( pThread-- > mThreads )
    {
        checkThreadOp(pthread_create(pThread,&mAttr,threadFunc,threadFuncArg),
                      "Thread creation failed: ");
    }
}

void ThreadPool::threadDeathMonitor( long newInterval,
                                     bool reportRestarts )
{
    if ( mNThreads )
    {
        long oldInterval;
        if ( true )
        {
            Locker ml(*this);

            oldInterval = mMonitorInterval;
            mMonitorInterval = newInterval;
            mReportRestarts = reportRestarts;
        }

        if ( !oldInterval && newInterval )
        {
            checkThreadOp(pthread_create(&mMonitorThread,0,threadFunc,this),
                          "Can't create ThreadPool monitor: ");
        }
        else if ( oldInterval && !newInterval )
        {
            mCondExit.signal();
            checkThreadOp(pthread_join(mMonitorThread,0),
                          "ThreadPool monitor join failed: ");
        }
    }
}

size_t ThreadPool::findThreadIndex( pthread_t handle )
{
    Locker ml(*this);
    pthread_t* pThread = mThreads + mNThreads;
    while ( pThread-- > mThreads )
    {
        if ( pthread_equal(handle,*pThread) )
            return pThread - mThreads;
    }
    return mNThreads;
}

void ThreadPool::shutdown()
{
    threadDeathMonitor(0);

    bool somebodyDied = false;
    pthread_t* pThread = mThreads + mNThreads;
    while ( pThread-- > mThreads )
    {
        void* retVal = 0;
        checkThreadOp(pthread_join(*pThread,&retVal),"Thread join failed: ");
        if ( reinterpret_cast<long>(retVal) )
        {
            std::cout << "Thread #" << (pThread-mThreads)
                        << " died with a return value of "
                        << reinterpret_cast<long>(retVal) << std::endl;
            somebodyDied = true;
        }
    }
    if ( somebodyDied )
        FatalErr("Quitting due to failure of child threads.");

    delete [] mThreads;
    mThreads = 0;
    mNThreads = 0;
}

void ThreadPool::monitorThreads()
{
    Locker ml(*this);
    while ( true )
    {
        ml.timedWait(mCondExit,mMonitorInterval);
        if ( !mMonitorInterval )
            break;

        pthread_t* pThread = mThreads + mNThreads;
        while ( pThread-- > mThreads )
        {
            if ( pthread_kill(*pThread,0) ==  ESRCH )
            {
                void* retVal;
                checkThreadOp(pthread_join(*pThread,&retVal),
                              "Can't join to dead thread: ");
                if ( mReportRestarts )
                    std::cout << "Restarting dead thread #" <<
                                 (pThread-mThreads) <<
                                 " which died with a return value of "
                                 << reinterpret_cast<long>(retVal) << std::endl;
                checkThreadOp(pthread_create(pThread,&mAttr,mThreadFunc,mThreadFuncArg),
                              "Can't create replacement thread: ");
            }
        }
    }
}

size_t ThreadPool::validateNThreads( size_t nThreads )
{
    size_t nProcs = processorsOnline();
    if ( nThreads > 3*nProcs )
        FatalErr("Trying to start a thread pool with " << nThreads <<
                 " threads, but there are only " << nProcs <<
                 " processors, which probably isn't going to work out well.");

    return nThreads;
}

void* ThreadPool::threadFunc( void* ptr )
{
    reinterpret_cast<ThreadPool*>(ptr)->monitorThreads();
    return 0;
}




/* -------------------------------------------------------------------------------
 * NaiveThreadPool
 *
 *  Associates an idle thread_ID with an item_ID.
 *
 * 2009-12 ribeiro@broadinstitute.org
 */

NaiveThreadPool::NaiveThreadPool(const int num_threads,
                                 const int num_items,
                                 pthread_mutex_t * mutex_lock) :
  num_threads_(num_threads),
  num_items_(num_items),
  mutex_lock_(mutex_lock)
{
  item_IDs_ = new int[num_threads];
  for (int i = 0; i < num_threads; i++)
    item_IDs_[i] = -1;

  thread_IDs_ = new int[num_items];
  for (int i = 0; i < num_items; i++)
    thread_IDs_[i] = -1;
}

NaiveThreadPool::~NaiveThreadPool()
{
  delete [] item_IDs_;
  delete [] thread_IDs_;
}

// Obtain a thread_ID and associate it with an item_ID
//    item_IDs_[thread_ID] gives the item_ID that thread_ID is processing
//    thread_IDs_[item_ID] gives the thread_ID that processed the item_ID
int NaiveThreadPool::AssignThread(const int item_ID)
{
  int thread_ID = 0;

  pthread_mutex_lock(mutex_lock_);

  // search list of thread_IDs for one that is available
  // (i.e., item_IDs_[thread_ID_] == -1)
  while (thread_ID < num_threads_ && item_IDs_[thread_ID] >= 0)
    thread_ID++;

  // if this assert fails something is wrong.
  // if we ask for a thread then there should be at least one available.
  ForceAssertLt(thread_ID, num_threads_);

  item_IDs_[thread_ID] = item_ID;
  thread_IDs_[item_ID] = thread_ID;
  pthread_mutex_unlock(mutex_lock_);

  return thread_ID;
}


void NaiveThreadPool::UnassignThread(const int item_ID)
{
  const int & thread_ID = thread_IDs_[item_ID];
  if (thread_ID != -1)
    item_IDs_[thread_ID] = -1;
  else
    std::cout << "trying to dissociate unassociated item_ID." << std::endl;
}








/* -----------------------------------------------------------------------------
 * ThreadBlockLocks
 *
 *  Multithreaded modules sometimes need to write to a shared data structure.
 *  You can mutex_lock everything down to write, but then the rest of the
 *  threads have to wait. The solution is to divide the data structure in blocks
 *  and lock only the block that is currently being written to, letting the
 *  other threads do something else. Also, data should be buffered so that you
 *  are not locking and unlocking all the time.
 *
 * 2009-12 ribeiro@broadinstitute.org
 */
ThreadBlockLocks::ThreadBlockLocks(const int num_blocks,
                                   const int num_threads,
                                   pthread_mutex_t * mutex_lock) :
  num_blocks_(num_blocks),
  num_threads_(num_threads),
  mutex_lock_(mutex_lock)
{
  thread_IDs_ = new int[num_blocks];
  for (int i = 0; i < num_blocks; i++)
    thread_IDs_[i] = -1;
}

ThreadBlockLocks::~ThreadBlockLocks()
{
  delete [] thread_IDs_;
}

// a debug function
void ThreadBlockLocks::OutputBlocked(const int block_ID, const int thread_ID)
{
  std::cout << "blocker("<< block_ID <<"): ";
  for (int i = 0; i != num_blocks_; i++) {
    const int & tID = thread_IDs_[i];
    if (tID < 0) std::cout << " -" << std::endl;
    else         std::cout << " " << tID << std::endl;
  }
}


// Locks block_ID, updating thread_IDs vec with the blocking thread_ID.
// returns the number of tries.
int ThreadBlockLocks::LockBlock(const int block_ID, const int thread_ID)
{
  int tries = 0;
  while (thread_IDs_[block_ID] != thread_ID) {

    pthread_mutex_lock(mutex_lock_);  // mutex_lock to access thread_IDs safely

    if (thread_IDs_[block_ID] < 0) // test for availability
      thread_IDs_[block_ID] = thread_ID;

    pthread_mutex_unlock(mutex_lock_);
    tries++;
  }
  return tries;
}

// Locks one block out of several possible in a set,
// updating thread_IDs vec with the blocking thread_ID.
// returns the block_ID.
int ThreadBlockLocks::LockSomeBlock(const std::set<int> & block_IDs,
                                    const int thread_ID)
{
  bool found = false;

  std::set<int>::iterator i_block_ID = block_IDs.begin();
  while (!found) {

    if (thread_IDs_[*i_block_ID] < 0) {   // block is not locked
      pthread_mutex_lock(mutex_lock_); // mutex_lock to access thread_IDs safely

      if (thread_IDs_[*i_block_ID] < 0) { // test for availability again
        thread_IDs_[*i_block_ID] = thread_ID;
        found = true;
      }

      pthread_mutex_unlock(mutex_lock_);
    }

    if (!found) {
      i_block_ID++;
      if (i_block_ID == block_IDs.end())
        i_block_ID = block_IDs.begin();
    }
  }
  return *i_block_ID;
}


// Unlocks block_ID.  Should only be done by the thread that locked the block.
// If so, there is no need to use mutex_locks.
void ThreadBlockLocks::UnlockBlock(const int block_ID, const int thread_ID)
{
  ForceAssertEq(thread_IDs_[block_ID], thread_ID);
  thread_IDs_[block_ID] = -1;
}

