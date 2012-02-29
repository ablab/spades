/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef BERNOULLI_H
#define BERNOULLI_H

#include "CoreTools.h"

// PartialBernoulliSum( n, k ): return sum_{i=0..k} choose(n,i).
//
// No attempt has been made to make this efficient or to pay attention to
// accuracy or overflow problems.

double PartialBernoulliSum( int n, int k );

// SurprisingTosses: Given a sequence of fair coin tosses, this looks at each 
// contiguous subsequence, and assigns it a p value, as follows: if it has length
// n and k heads in it, the p value is PartialBernoulliSum( n, Min( k, n-k ) )/2^n.
// The return value is the minimum of all these p values.
//
// No attempt has been made to make this efficient or to pay attention to
// accuracy or overflow problems.

double SurprisingTosses( const vec<Bool>& s, int max_seq = -1 );

// BinomialSum.  Compute sum_{i=0}^k (n choose i) p^i (1-p)^(n-i).
// Note that class T needs to be a class that can hold rational numbers,
// e.g. float or double.  There is a special version for T = double that is
// implemented slightly differently.

double BinomialSum( int n, int k, double p );

template<class T> T BinomialSum( int n, int k, T p )
{    ForceAssertGe( n, 1 );
     ForceAssertGe( k, 0 );
     ForceAssertLe( k, n );
     T sum = 0, choose = 1, product = IPow( T(1-p), n );
     for ( int i = 0; i <= k; i++ )
     {    sum += choose * product;
          choose *= T(n - i);
          choose /= T(i + 1);
          product *= p / (1 - p);    }
     return sum;    }

#endif
