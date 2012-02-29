///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file Primes.h
 * \author tsharpe
 * \date Feb 4, 2011
 *
 * \brief A class to group functions dealing with powers of 2.
 */
#ifndef PRIMES_H_
#define PRIMES_H_

#include "system/Assert.h"

/// Some functions relating to powers of 2
class PowerOf2
{
public:
    /// number of leading 0 bits
    static int nlz( unsigned long val )
    { return val ? __builtin_clzl(val) : 64; }

    /// return floor(log2(val)).  floorLg2(0) is -1.
    static int floorLg2( unsigned long val )
    { return 63-nlz(val); }

    /// return ceil(log2(val)). ceilLg2(0) is -1.
    static int ceilLg2( unsigned long val )
    { return val ? 64-nlz(val-1) : -1; }

    /// return greatest power of 2 <= val.  floor2(0) is 0.
    static unsigned long floor2( unsigned long val )
    { return val ? 1ul << floorLg2(val) : 0ul; }

    /// return least power of 2 >= val.  ceil2(0) is 0.
    static unsigned long ceil2( unsigned long val )
    { return val > 1 ? 1ul << ceilLg2(val) : 0ul; }

    /// Return the largest prime number <= a specified power of 2.
    /// Good up to 2^64.  getNearbyPrime(0) is 1.
    static unsigned long getNearbyPrime( unsigned powerOf2 );
};

#endif /* PRIMES_H_ */
