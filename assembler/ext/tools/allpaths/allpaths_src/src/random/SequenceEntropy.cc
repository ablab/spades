/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2010) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"

vec<int> DEcount(16);

vec<double> Logs( int n )
{    vec<double> L(n);
     for ( int i = 1; i < n; i++ )
          L[i] = log( double(i) );
     return L;    }

vec<double> ilog( Logs(200) );

double DinukeEntropy( const basevector& b )
{    for ( int i = 0; i < 16; i++ )
          DEcount[i] = 0;
     for ( int i = 0; i < b.isize( ) - 1; i++ )
          ++DEcount[ ( 4 * b[i] ) + b[i+1] ];
     double sum = 0.0;
     for ( int i = 0; i < 16; i++ )
     {    if ( DEcount[i] == 0 ) continue;
          double p = double( DEcount[i] ) / double( b.size( ) - 1 );
          double logp;
          if ( DEcount[i] < ilog.isize( ) && b.isize( ) - 1 < ilog.isize( ) )
               logp = ilog[ DEcount[i] ] - ilog[ b.isize( ) - 1 ];
          else logp = log(p);
          sum -= p * logp;    }
     return sum;    }
