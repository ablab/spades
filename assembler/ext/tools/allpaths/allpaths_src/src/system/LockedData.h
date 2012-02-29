///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file LockedData.h
 * \author tsharpe
 * \date Jul 17, 2009
 *
 * \brief
 */

#ifndef SYSTEM_LOCKEDDATA_H_
#define SYSTEM_LOCKEDDATA_H_

#include "system/Assert.h"
#include <cerrno>
#include <pthread.h>

/// A class that contains a mutex.
/// The intent is that you'll extend the class to include the
/// data controlled by the mutex.
class LockedData
{
public:
    LockedData()
    { checkThreadOp(pthread_mutex_init(&mMutex, 0),
                    "Mutex initialization failed."); }

    ~LockedData()
    { checkThreadOp(pthread_mutex_destroy(&mMutex),
                    "Mutex destruction failed."); }

    void lock()
    { checkThreadOp(pthread_mutex_lock(&mMutex), "Mutex locking failed."); }

    void unlock()
    { checkThreadOp(pthread_mutex_unlock(&mMutex), "Mutex unlocking failed."); }

    /// prints an error message and exits if errNo is not 0
    static void checkThreadOp( int errNo, char const* msg )
    { if ( errNo ) die(errNo,msg); }

private:
    LockedData( LockedData const& ); // not implemented -- no copying
    LockedData& operator=( LockedData const& ); // not implemented -- no copying

    static void die( int errNo, char const* msg );

    pthread_mutex_t mMutex;
    friend class Condition;
};

/// A class to encapsulate condition variables.
class Condition
{
public:
    Condition( LockedData& lock ) : mLock(lock)
    { LockedData::checkThreadOp(pthread_cond_init(&mCV,0),
                                "Unable to initialize condition variable: "); }
    ~Condition()
    { LockedData::checkThreadOp(pthread_cond_destroy(&mCV),
                                "Unable to destroy condition variable: "); }

    void signal()
    { LockedData::checkThreadOp(pthread_cond_signal(&mCV),
                                "Unable to signal condition: "); }

    void broadcast()
    { LockedData::checkThreadOp(pthread_cond_broadcast(&mCV),
                                "Unable to broadcast condition: "); }

    /// Wait for condition.
    /// You must have the associated lock already locked. So consider calling
    /// this via a Locker, instead.
    void wait()
    { LockedData::checkThreadOp(pthread_cond_wait(&mCV,&mLock.mMutex),
                                "Unable to wait for condition: "); }

    /// Wait for condition.
    /// You must have the associated lock already locked.  Consider calling this
    /// from a Locker, rather than directly.
    /// Returns true if we timed out.
    bool timedWait( long nSecs )
    { timespec absTime; absTime.tv_sec = time(0) + nSecs; absTime.tv_nsec = 0L;
      int err = pthread_cond_timedwait(&mCV, &mLock.mMutex, &absTime);
      if ( err != ETIMEDOUT )
         LockedData::checkThreadOp(err,"Unable to wait for condition: ");
      return err == ETIMEDOUT; }

private:
    Condition( Condition const& ); // unimplemented -- no copying
    Condition& operator=( Condition const& ); // unimplemented -- no copying

    LockedData& mLock;
    pthread_cond_t mCV;

    friend class Locker;
};

/// A class to provide thread-safe access to data via a mutex.
/// The intent is that you'll extend this class with accessors for the
/// controlled data which you'll trivially forward to the LockedData
/// class.
class Locker
{
public:
    Locker( LockedData& lock ) : mLock(lock) { mLock.lock(); }
    ~Locker() { mLock.unlock(); }

    void wait( Condition& condition )
    { AssertEq(&mLock,&condition.mLock);
      condition.wait(); }

    bool timedWait( Condition& condition, long nSecs )
    { AssertEq(&mLock,&condition.mLock);
      return condition.timedWait(nSecs); }

protected:
    LockedData& mLock;

private:
    Locker( Locker const& ); // unimplemented -- no copying
    Locker& operator=( Locker const& ); // unimplemented -- no copying
};

#endif /* SYSTEM_LOCKEDDATA_H_ */
