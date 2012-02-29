///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file RNGen.h
 * \author tsharpe
 * \date Nov 17, 2011
 *
 * \brief
 */
#ifndef RANDOM_RNG_H_
#define RANDOM_RNG_H_

#include "system/SpinLockedData.h"
#include <algorithm>
#include <climits>

/// Random number generator.
/// Now that it's a class, you can have multiple independent sequences, and you
/// can save and restore state simply by copying RNGen objects.
/// The static methods random and srandom are thread-safe (via a spin lock):
/// you'd be better off having a separate RNGen for each thread, but sometimes
/// that's hard to arrange.
class RNGen
{
    typedef unsigned long state_t;
public:
    RNGen() { *this = gDflt; }
    RNGen( unsigned seedVal ) { seed(seedVal); }
    RNGen( RNGen const& that ) { *this = that; }

    // compiler-supplied destructor is OK

    RNGen& operator=( RNGen const& that )
    { std::copy(that.mState,that.mState+STATE_SIZE,mState);
      mpFront = mState + (that.mpFront-that.mState);
      mpRear = mState + (that.mpRear-that.mState);
      return *this; }

    long next()
    { unsigned result = (*mpFront += *mpRear);
      if ( ++mpFront >= mState+STATE_SIZE )
      { mpFront = mState; ++mpRear; }
      else if ( ++mpRear >= mState+STATE_SIZE )
        mpRear = mState;
      return result >> 1; }

    void seed( unsigned seedVal )
    { state_t last = seedVal;
      mState[0] = last;
      state_t* ppp = mState;
      while ( ++ppp < mState+STATE_SIZE )
        *ppp = last = (last*1103515245 + 12345);
      mpFront = mState + 3;
      mpRear = mState;
      int nnn = 10*31;
      while ( nnn-- )
        next(); }

    static long random()
    { SpinLocker locker(gLock); return gSystem.next(); }

    static void srandom( unsigned seedVal )
    { SpinLocker locker(gLock); gSystem.seed(seedVal); }

    static long const RNGEN_RAND_MAX = INT_MAX;

private:
    static RNGen gDflt;
    static RNGen gSystem;
    static SpinLockedData gLock;

    static unsigned long const STATE_SIZE = 31;

    unsigned long mState[STATE_SIZE];
    unsigned long* mpFront;
    unsigned long* mpRear;
};

// Adapted from code that included this notice:
/*
 * Copyright (c) 1983 Regents of the University of California.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms are permitted
 * provided that the above copyright notice and this paragraph are
 * duplicated in all such forms and that any documentation,
 * advertising materials, and other materials related to such
 * distribution and use acknowledge that the software was developed
 * by the University of California, Berkeley.  The name of the
 * University may not be used to endorse or promote products derived
 * from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 */


#endif /* RANDOM_RNG_H_ */
