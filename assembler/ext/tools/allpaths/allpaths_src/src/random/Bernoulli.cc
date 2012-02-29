/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "math/Functions.h"
#include "random/Bernoulli.h"

double PartialBernoulliSum( int n, int k )
{    ForceAssertGe( n, 1 );
     ForceAssertGe( k, 0 );
     ForceAssertLe( k, n );
     double sum = 0;
     double choose = 1.0;
     for ( int i = 0; i <= k; i++ )
     {    sum += choose;
          choose *= double(n - i);
          choose /= double(i + 1);    }
     return sum;    }

double SurprisingTosses( const vec<Bool>& s, int max_seq )
{    double minp = 1.0;
     for ( int start = 0; start < s.isize( ); start++ )
     {    int top = s.size( );
          if ( max_seq >= 0 ) top = Min( top, start + max_seq );
          for ( int stop = start+1; stop < top; stop++ )
          {    int n = stop - start + 1, k = 0;
               for ( int j = start; j <= stop; j++ )
                    if ( s[j] ) ++k;
               if ( n-k < k ) k = n-k;
               double ppow = pow( 0.5, n );
               minp = Min( minp, PartialBernoulliSum( n, k ) * ppow );    }    }
     return minp;    }

double BinomialSum( int n, int k, double p )
{    ForceAssertGe( n, 1 );
     ForceAssertGe( k, 0 );
     ForceAssertLe( k, n );
     long double sum = 0, choose = 1.0, product = pow( 1.0 - p, double(n) );
     for ( int i = 0; i <= k; i++ ) 
     {    sum += choose * product;
          choose *= double(n - i);
          choose /= double(i + 1);
          product *= p / (1.0 - p);    }
     return sum;    }
