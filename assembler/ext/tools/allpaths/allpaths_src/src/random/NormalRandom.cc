/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "random/NormalRandom.h"
#include "random/Random.h"

// class NormalRandom
//
// Provides a source of random numbers with a normal distribution.
//

double
NormalRandom::value() const
{
  // method for getting the next value of the normal random variable  

  const long two_pwr_30 = 1073741824;
  const double two_pwr_31 = 2147483648.0;
  const double v_denom = 117732250.0;
  
  bool done = false;
  double answer = 0.0;
  while (!done)
    {
      // Note: randomx() returns a long value uniformly distributed in [0,2^31-1]
      long x = randomx();
      long y = randomx();

      // u is a random variable uniformly distributed in (0,1]
      // v is a random variable uniformly distributed in (-V,+V),
      // where V is (a slight over-esitmation of)  2*sqrt(ln(2^30))
      double u = static_cast<double>(x) / two_pwr_31;
      double v = static_cast<double>(y - two_pwr_30) / v_denom;

      // avoid any chance that u may be "zero"
      if (u < 1e-8)
	continue;

      // call NormalDeviate(...) to get a random normal variable in answer;
      // redo everything if this call fails
      if (NormalDeviate(u, v, answer))
	done = true;
    }

  return mean_ + stddev_*answer;
}

double FastNormal( ) 
{    static vec<double> normals;
     if ( normals.empty( ) )
     {    normals.resize(1000000);
          NormalRandom r( 0, 1 );
          for ( int i = 0; i < normals.isize( ); i++ )
               normals[i] = r.value( );    }
     static int count(0); 
     if ( count == normals.isize( ) ) count = 0;
     return normals[count++];    }
