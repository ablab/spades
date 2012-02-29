///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file WorklistUtils.h
 * \author tsharpe
 * \date Feb 4, 2009
 *
 * \brief A couple of utility classes for the Worklist.
 */
#ifndef SYSTEM_WORKLISTUTILS_H_
#define SYSTEM_WORKLISTUTILS_H_


#include "system/LockedData.h"
#include "system/SpinLockedData.h"
#include <cstddef>
#include <iostream>
#include <pthread.h>
#include <set>

class WorklistParameterizer
{
public:
    WorklistParameterizer( size_t nUnits, size_t bytesPerUnit,
                           size_t minBatchSize, unsigned nThreads );

    // compiler-supplied copying and destructor are OK

    size_t getBatchSize() const { return mBatchSize; }
    size_t getNBatches() const { return mNBatches; }
    unsigned getNThreads() const { return mNThreads; }

private:
    size_t mBatchSize;
    size_t mNBatches;
    unsigned mNThreads;
};

class Dotter : private SpinLockedData
{
public:
    Dotter( size_t nBatches = 0 ) : mNBatches(nBatches), mBatch(0) {}

    // compiler-supplied copying and destructor are OK

    size_t getNBatches() const { return mNBatches; }
    void setNBatches( size_t nBatches ) { mNBatches = nBatches; mBatch = 0; }

    void batchDone()
    { SpinLocker locker(*this);
      if ( mNBatches <= 100 || 100*mBatch/mNBatches < 100*(mBatch+1)/mNBatches )
          std::cout << '.' << std::flush;
      if ( (mBatch += 1) == mNBatches ) std::cout << std::endl; }

private:
    size_t mNBatches;
    size_t mBatch;
};

/// A class to track the state of the queue.
class QueueState : public LockedData
{
public:
    QueueState()
    : mSize(0), mDone(false), mCondNotEmpty(*this), mCondEmpty(*this) {}

    // The entire interface, except for the constructor and destructor, is private because
    // the state should only be checked when the mutex is locked. The QueueStateManipulator,
    // which is a friend, locks the mutex in its constructor.
private:
    void setDone() { mDone = true; }
    bool isDone() { return mDone; }

    size_t incSize() { return ++mSize; }
    size_t decSize() { return --mSize; }
    size_t getSize() { return mSize; }
    void setSize( size_t sz ) { mSize = sz; }

    QueueState( QueueState const& ); // unimplemented -- no copying
    QueueState& operator=( QueueState const& ); // unimplemented -- no copying

    size_t mSize;            // current size of the work queue
    bool mDone;              // nothing more will be added to the queue
    Condition mCondNotEmpty; // predicate "queue not empty"
    Condition mCondEmpty;    // predicate "queue is empty"
    friend class QueueStateManipulator;
};

/// A class to change the queue state.
/// Constructor locks the QueueState mutex.
/// Destructor unlocks and does any necessary condition signaling.
class QueueStateManipulator
{
public:
    QueueStateManipulator( QueueState& qs )
    : mQS(qs)
    {
        mQS.lock();
        mSize = mQS.getSize();
        mDone = mQS.isDone();
    }

    ~QueueStateManipulator();

    void setDone() { mQS.setDone(); }
    bool isDone() { return mQS.isDone(); }

    size_t incSize() { return mQS.incSize(); }
    size_t decSize() { return mQS.decSize(); }
    size_t getSize() { return mQS.getSize(); }
    void setSize( size_t sz ) { mQS.setSize(sz); }

    void waitForEmpty()
    { while ( getSize() ) mQS.mCondEmpty.wait(); }

    void waitForWork()
    { while ( !getSize() && !isDone() ) mQS.mCondNotEmpty.wait(); }

private:
    QueueStateManipulator( QueueStateManipulator const& ); // unimplemented -- no copying
    QueueStateManipulator& operator=( QueueStateManipulator const& ); // unimplemented -- no copying

    QueueState& mQS;
    size_t mSize; // size of queue at constructor-time
    bool mDone;   // state of done flag at constructor-time
};

/// Class to manage a pool of threads
class ThreadPool : LockedData
{
public:
    ThreadPool( size_t nThreads, size_t threadStackSize,
                void* (*threadFunc)(void*), void* threadFuncArg );

    ~ThreadPool()
    { shutdown();
      checkThreadOp(pthread_attr_destroy(&mAttr),
                    "ThreadAttr destruction failed: "); }


    size_t getNThreads() { return mNThreads; }

    void threadDeathMonitor( long howOftenToCheckInSecs,
                             bool reportRestarts = true );
    size_t findThreadIndex( pthread_t );

    void shutdown();

private:
    void monitorThreads();
    static size_t validateNThreads( size_t nThreads );
    static void* threadFunc( void* ptr );

    void* (*mThreadFunc)(void*);
    void* mThreadFuncArg;
    size_t mNThreads;
    pthread_t* mThreads;
    pthread_t mMonitorThread;
    pthread_attr_t mAttr;
    Condition mCondExit;    // predicate "queue is empty"
    long mMonitorInterval;
    bool mReportRestarts;
    static size_t const DEFAULT_THREAD_STACK_SIZE = 8ul*1024ul*1024ul;
};


/* -------------------------------------------------------------------------------
 * NaiveThreadPool
 * 
 *  Associates an idle thread_ID with an item_ID.
 *
 * 2009-12 ribeiro@broadinstitute.org 
 */

class NaiveThreadPool
{
private:
  int num_threads_;
  int num_items_;
  pthread_mutex_t * mutex_lock_;

  int * item_IDs_;
  int * thread_IDs_;

public:
  NaiveThreadPool(const int num_threads, 
                  const int num_items, 
                  pthread_mutex_t * mutex_lock);

  ~NaiveThreadPool(); 

  // Obtain a thread_ID and associate it with an item_ID
  //    item_IDs_[thread_ID] gives the item_ID that thread_ID is processing
  //    thread_IDs_[item_ID] gives the thread_ID that processed the item_ID
  int AssignThread(const int item_ID);
  void UnassignThread(const int item_ID);

  inline int GetNumThreads() const { return num_threads_; }
  inline int GetNumItems() const { return num_items_; }
};




/* -------------------------------------------------------------------------------
 * ThreadBlockLocks
 *
 *  Multithreaded modules sometimes need to write to a shared data structure.
 *  You can mutex_lock everything down to write, but then the rest of the threads have to wait.
 *  The solution is to divide the data structure in blocks and lock only the 
 *  block that is currently being written to, letting the other threads do something else.    
 *  Also, data should be buffered so that you are not locking and unlocking all the time.
 *
 * 2009-12 ribeiro@broadinstitute.org 
 */
class ThreadBlockLocks
{
private:
  int num_blocks_;
  int num_threads_;
  pthread_mutex_t * mutex_lock_;

  int * thread_IDs_;
  
public:
  ThreadBlockLocks(const int num_blocks, 
                   const int num_threads, 
                   pthread_mutex_t * mutex_lock);
  
  ~ThreadBlockLocks();

  inline int GetNumBlocks() const { return num_blocks_; }
  inline int GetNumThreads() const { return num_threads_; }
  
  // a debug function
  void OutputBlocked(const int block_ID, const int thread_ID);
  
  // Locks block_ID, updating thread_IDs vec with the blocking thread_ID.
  // returns the number of tries. 
  int LockBlock(const int block_ID, const int thread_ID);

  // Locks one block out of several possible in a set, 
  // updating thread_IDs vec with the blocking thread_ID.
  // returns the block_ID. 
  int LockSomeBlock(const std::set<int> & block_IDs, const int thread_ID);

  // Unlocks block_ID.  Should only be done by the thread that locked the block.  
  // If so, there is no need to use mutex_locks.
  void UnlockBlock(const int block_ID, const int thread_ID);
};













#endif /* SYSTEM_WORKLISTUTILS_H_ */
