///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// AddUnibaseJunctions.  Add (K+1)-mers to unibases, each of which represents a 
// junction between unibases.  Note that the output "unibases" are no longer 
// true unibases.

#include "Basevector.h"
#include "MainTools.h"
#include "paths/GetNexts.h"

int main( int argc, char** argv ) 
{
  
     RunTime();

     BeginCommandArguments;
     CommandArgument_String( PRE );
     CommandArgument_String( DATA );
     CommandArgument_String( RUN );
     CommandArgument_Int( K );
     CommandArgument_String_OrDefault( UNIBASES_IN, "all_reads" );
     CommandArgument_String_OrDefault( UNIBASES_OUT, UNIBASES_IN + ".junctionized" );
     EndCommandArguments;

     // Define directories.
  
     String data_dir = PRE + "/" + DATA;
     String run_dir = data_dir + "/" + RUN;

     // Load unibases.

     String kK = ".k" + ToString(K);
     vecbasevector unibases( run_dir + "/" + UNIBASES_IN + ".unibases" + kK );

     // Add junctions.

     vec< vec<int> > nexts;
     GetNexts( K, unibases, nexts );
     size_t nuni = unibases.size( ), add = 0;
     for ( size_t u = 0; u < nuni; u++ )
          add += nexts[u].size( );
     unibases.reserve( nuni + add );
     basevector b;
     for ( size_t u = 0; u < nuni; u++ )
     {    for ( int j = 0; j < nexts[u].isize( ); j++ )
          {    b = unibases[u];
               b.resize( b.size( ) + 1 );
               b.Set( b.size( ) - 1, unibases[ nexts[u][j] ][K-1] );
               unibases.push_back(b);    }    }

     // Write new unibases.

     unibases.WriteAll( run_dir + "/" + UNIBASES_OUT + ".unibases" + kK );    }
