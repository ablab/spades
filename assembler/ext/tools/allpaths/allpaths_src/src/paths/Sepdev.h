///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SEPDEV_H
#define SEPDEV_H

#include "CoreTools.h"

// Class template: Tsepdev
//
// A sepdev defines a separation between two unspecified things, along with a 
// deviation value for it.
template<typename T>
class Tsepdev {

     public:

     typedef T value_type;

     Tsepdev( ) { }
     Tsepdev( T sep, T dev ) : sep_(sep), dev_(dev) { }
     template<typename U>
     explicit Tsepdev( U sep, U dev ) : sep_(sep), dev_(dev) { }
     // specialization later for T=int, U=float/double

     T Sep( ) const { return sep_; }
     T Dev( ) const { return dev_; }

     // things that make it easier to modify them:
     void Flip() { sep_ = -sep_; }
     void AddToSep( T t ) { sep_ += t; }

     private:

     T sep_;
     T dev_;

};

typedef Tsepdev<int> sepdev;
typedef Tsepdev<double> fsepdev;

// Templatized constructor helpers to convert to int intelligently.
template<> template<> 
inline Tsepdev<int>::Tsepdev( float sep, float dev )
  : sep_(int(round(sep))), dev_(int(ceil(dev))) { }
template<> template<> 
inline Tsepdev<int>::Tsepdev( double sep, double dev )
  : sep_(int(round(sep))), dev_(int(ceil(dev))) { }

void BinaryWrite( int fd, const sepdev& sd );
void BinaryRead( int fd, sepdev& sd );
void BinaryWrite( int fd, const fsepdev& fsd );
void BinaryRead( int fd, fsepdev& fsd );

#endif

