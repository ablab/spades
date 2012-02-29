///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <vector>

#include "math/Arith.h"

int Ilog10( double x )
{    Bool first_call(True);
     const int max_power = 100;
     static vector<double> ten_to( max_power + 1 );
     if (first_call)
     {    ten_to[0] = 1.0;
          for ( int i = 1; i <= max_power; i++ )
               ten_to[i] = ten_to[i-1] * 10.0;
          first_call = False;    }
     int answer = int( floor( log10(x) ) );
     if ( answer >= 0 && answer <= max_power )
     {    if ( ten_to[answer] > x ) --answer;
          else if ( answer < max_power && ten_to[answer+1] <= x ) ++answer;    }
     return answer;    }

