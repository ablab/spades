///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "random/Random.h"
#include "reporting/ScaffoldAccuracyCore.h"

void SelectTestSequences( int& SAMPLE, const int D, const vecbasevector& assembly, 
     const vecbitvector& assembly_amb, vecbasevector& query,
     vec<int>& t1s, vec<int>& t2s, vec<int>& start1s, vec<int>& start2s )
{
     query.clear( ), t1s.clear( ), t2s.clear( ), start1s.clear( ), start2s.clear( );
     size_t A = assembly.sumSizes();
     const int n = 100;
     if ( A-D-n <= 0 )
     {    cout << "Can't find SAMPLE=" << SAMPLE << " test cases." << endl;
          SAMPLE = 0;
          cout << "Resetting SAMPLE to 0." << endl;    }
     else
     {    int tries = 0;
          for ( int i = 0; i < SAMPLE; i++ )
          {    int start1 = randomx( ) % (A-D-n);
               int start2 = start1 + D;
               size_t t1, t2;
               for ( t1 = 0; t1 < assembly.size( ); t1++ )
               {    if ( start1 < assembly[t1].isize( ) ) break;
                    start1 -= assembly[t1].size( );    }
	       if ( t1 == assembly.size() ) continue;
	       
               for ( t2 = 0; t2 < assembly.size( ); t2++ )
               {    if ( start2 < assembly[t2].isize( ) ) break;
                    start2 -= assembly[t2].size( );    }
	       if ( t2 == assembly.size() ) continue;
	       
               if ( t1 != t2 || start2 > assembly[t2].isize( ) - n )
               {    --i;
                    if ( ++tries > 10*SAMPLE )
		      {    cout << "Can't find SAMPLE=" << SAMPLE << " test cases." << endl;
                         SAMPLE = i+1;
                         cout << "Resetting SAMPLE to " << SAMPLE << "." << endl;
                         break;    }
                    continue;    }
               Bool ambig = False;
               for ( int j = 0; j < n; j++ )
               {    if ( assembly_amb[t1][start1+j] ) ambig = True;
                    if ( assembly_amb[t2][start2+j] ) ambig = True;    }
               if (ambig)
               {    --i;
                    continue;    }
               t1s.push_back(t1), t2s.push_back(t2);
               start1s.push_back(start1), start2s.push_back(start2);
               basevector e1, e2;
               e1.SetToSubOf( assembly[t1], start1, n );
               e2.SetToSubOf( assembly[t2], start2, n );
               query.push_back_reserve(e1); 
               query.push_back_reserve(e2);    }    }    }
