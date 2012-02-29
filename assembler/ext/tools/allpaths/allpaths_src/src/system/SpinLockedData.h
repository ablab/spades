///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file SpinLockedData.h
 * \author tsharpe
 * \date Nov 16, 2011
 *
 * \brief
 */
#ifndef SYSTEM_SPINLOCKEDDATA_H_
#define SYSTEM_SPINLOCKEDDATA_H_

#include "system/Assert.h"

/// A spin-lock.
class SpinLockedData
{
public:
    SpinLockedData() : mLockByte(0) {}

    ~SpinLockedData() { AssertEq(mLockByte,0); }

    void lock()
    { while ( __sync_lock_test_and_set(&mLockByte,1) )
        while ( mLockByte )
        {} }

    void unlock()
    { AssertEq(mLockByte,1); __sync_lock_release(&mLockByte); }

private:
    SpinLockedData( SpinLockedData const& ); // not implemented -- no copying
    SpinLockedData& operator=( SpinLockedData const& ); // not implemented -- no copying

    char volatile mLockByte;
};

/// Something that operates a spin-lock, and never forgets to unlock it.
class SpinLocker
{
public:
    SpinLocker( SpinLockedData& lock ) : mLock(lock) { mLock.lock(); }
    ~SpinLocker() { mLock.unlock(); }

private:
    SpinLocker( SpinLocker const& ); // unimplemented -- no copying
    SpinLocker& operator=( SpinLocker const& ); // unimplemented -- no copying

    SpinLockedData& mLock;
};

#endif /* SYSTEM_SPINLOCKEDDATA_H_ */
