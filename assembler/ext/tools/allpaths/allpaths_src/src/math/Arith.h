///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Define a class Float whose sole purpose is to hold a float and facilitate
// platform-independent reproducibility of computations on floats.  
//
// Also define Ilog10.

#ifndef ARITH_H
#define ARITH_H

#include "system/Types.h"
#include "system/System.h"
#include <math.h>
#include <iostream>

class Float {
     public:
     Float( ) { }
     Float( float x ) { value_ = x; }
     float Value( ) const { return value_; }
     operator float( ) { return value_; }
     void Set( float x ) { value_ = x; }

     void operator+=( Float x )
     {    
          #ifdef i386
               Set( float( double(Value( )) + double(x.Value( )) ) );
          #else
               Set( Value( ) + x.Value( ) );
          #endif
               }

     void operator-=( Float x )
     {    
          #ifdef i386
               Set( float( double(Value( )) - double(x.Value( )) ) );
          #else
               Set( Value( ) - x.Value( ) );
          #endif
               }

     void operator*=( Float x )
     {    
          #ifdef i386
               Set( float( double(Value( )) * double(x.Value( )) ) );
          #else
               Set( Value( ) * x.Value( ) );
          #endif
               }

     void operator/=( Float x )
     {    
          #ifdef i386
               Set( float( double(Value( )) / double(x.Value( )) ) );
          #else
               Set( Value( ) / x.Value( ) );
          #endif
               }

     friend Bool operator<( Float x, Float y )
     {    
          #ifdef i386
               return (double) x.Value( ) < (double) y.Value( );
          #else
               return x.Value( ) < y.Value( );
          #endif
               }

     friend Bool operator>( Float x, Float y )
     {    return y < x;    }

     friend Bool operator>=( Float x, Float y )
     {    return !( x < y );    }

     friend Bool operator<=( Float x, Float y )
     {    return !( y < x );    }

     friend void BinaryWrite( int fd, const Float& f ) {
       WriteBytes( fd, &f.value_, sizeof(f.value_) );
     }
     
     friend void BinaryRead( int fd, Float& f ) {
       ReadBytes( fd, &f.value_, sizeof(f.value_) );
     }

     private:
     float value_;
};

inline Float operator+( Float x, Float y )
{    
     #ifdef i386
	  return float(double(x.Value( ))+double(y.Value( )));
     #else
          return x.Value( )+y.Value( );
     #endif
          }

inline Float operator-( Float x, Float y )
{    
     #ifdef i386
	  return float(double(x.Value( ))-double(y.Value( )));
     #else
          return x.Value( )-y.Value( );
     #endif
          }

inline Float operator*( Float x, Float y )
{    
     #ifdef i386
	  return float(double(x.Value( ))*double(y.Value( )));
     #else
          return x.Value( )*y.Value( );
     #endif
          }

inline Float operator/( Float x, Float y )
{    
     #ifdef i386
	  return float(double(x.Value( ))/double(y.Value( )));
     #else
          return x.Value( )/y.Value( );
     #endif
          }

inline Float Pow( Float x, Float y )
{    return Float(static_cast<float>( pow( x.Value( ), y.Value( ) ) ) );    }

inline Float Pow10( Float y )
{    return Pow( Float(10.0), y );    }

// Ilog10(x) returns the smallest integer power of 10 which is <= x.  It makes
// a feeble attempt to guarantee correctness.

int Ilog10( double x );

inline
ostream& operator<< ( ostream& out, const Float& aFloat ) {
  return out << aFloat.Value();
}

inline
istream& operator>> ( istream& in, Float& aFloat ) {
  float f;
  in >> f;
  aFloat.Set( f );
  return in;
}

#endif
