///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Class and functions related to one-dimensional array

#ifndef ARRAY_H
#define ARRAY_H

#include "CoreTools.h"
#include "Vec.h"

// Smooth an array of double values using an Gaussian. 
inline void SmoothArrayGaussian ( vec<double> & distr, int delta )
{
     // prepare the kernal function for smoothing
     vec<double> kernal( 8 * delta + 1 );
     int center = 4 * delta; // center position of the gaussian
     for ( int i = 0; i < kernal.isize(); i++ ) {
          int x = i - center;
          kernal[i] = exp( - x * x / double( 2 * delta * delta ) );
     }
     double total = Sum( kernal );
     for ( int i = 0; i < kernal.isize(); i++ ) 
          kernal[i] /= total;
     // smoothed array
     vec<double> distr_s( distr.size(), 0 );
     for ( int i = 0; i < distr.isize(); i++ ) 
          for ( int k = 0; k < kernal.isize(); k++ ) {
               int x = i + k - center;
               if ( x < 0 || x >= distr.isize() ) continue;
               distr_s[x] += distr[i] * kernal[k];
          }
     swap(distr, distr_s);
}

// Smooth an array of double values using an Gaussian. 
// Here the array is sparse, and is represented using std::map<int, double>
inline void  SmoothArrayGaussian ( map<int, double> & distr, int delta )
{
     // prepare the kernal function for smoothing
     vec<double> kernal( 8 * delta + 1 );
     int center = 4 * delta; // center position of the gaussian
     for ( int i = 0; i < kernal.isize(); i++ ) {
          int x = i - center;
          kernal[i] = exp( - x * x / double( 2 * delta * delta ) );
     }
     double total = Sum( kernal );
     for ( int i = 0; i < kernal.isize(); i++ ) 
          kernal[i] /= total;
     // smoothed array
     map<int, double> distr_s;
     for ( map<int, double>::iterator it = distr.begin(); 
               it != distr.end(); it++ )
          for ( int k = 0; k < kernal.isize(); k++ ) {
               int x = it->first  + k - center;
               distr_s[x] += it->second * kernal[k];
          }
     swap(distr, distr_s);
}

// Smoothing with an output
template <typename T>
void SmoothArrayGaussian ( const T& distr, int delta,
         T & distr_s )
{
     distr_s = distr;
     SmoothArrayGaussian( distr_s, delta );
}
#endif
