///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// GetNexts: determine which unibase can follow which, assuming that a K-1 base
// overlap is good enough.  Can be applied to things other than unibases.

#ifndef GET_NEXTS
#define GET_NEXTS

#include "Basevector.h"
#include "CoreTools.h"
#include "kmers/KmerRecord.h"

// T is intended to be either a vecbasevector or vec<basevector>

template<int K, class T> void GetNexts( const T& unibases, vec< vec<int> >& nexts )
{    int nuni = unibases.size( );
     nexts.clear_and_resize(nuni);
     {    vec< kmer_record<K-1,1> > unistarts, unistops;
          unistarts.reserve(nuni), unistops.reserve(nuni);
          for ( size_t i = 0; i < unibases.size( ); i++ )
          {    basevector b;
               b.SetToSubOf( unibases[i], 0, K-1 );
               kmer_record<K-1,1> rec;
               rec.Set( b, i, 0 );
               unistarts.push_back(rec);
               int n = unibases[i].size( );
               b.SetToSubOf( unibases[i], n - (K-1), K-1 );
               rec.Set( b, i, n - (K-1) );
               unistops.push_back(rec);    }
          Sort(unistarts), Sort(unistops);
          int ulast = 0;
          for ( int i = 0; i < unistops.isize( ); i++ )
          {    int j;
               Bool eq = False;
               for ( j = ulast; j < unistarts.isize( ); j++ )
               {    if ( unistarts[j] > unistops[i] ) break;
                    if ( unistarts[j].EqualKmers( unistops[i] ) )
                    {    eq = True;
                         break;    }    }
               if ( j == unistarts.isize( ) ) break;
               if (eq)
               {    for ( int z = j; z < unistarts.isize( ); z++ )
                    {    if ( !unistarts[z].EqualKmers( unistops[i] ) ) break;
                         nexts[ unistops[i].GetId( ) ].
                              push_back( unistarts[z].GetId( ) );    }    }
               ulast = j;    }    }    }



// Wrapper function to handle templatization.

template<class T> void GetNexts( const int K, const T & unibases,
	       vec< vec<int> > & nexts )
{
  switch ( K ) {
  case 16: GetNexts<16>( unibases, nexts ); break;
  case 20: GetNexts<20>( unibases, nexts ); break;
  case 40: GetNexts<40>( unibases, nexts ); break;
  case 60: GetNexts<60>( unibases, nexts ); break;
  case 64: GetNexts<64>( unibases, nexts ); break;
  case 80: GetNexts<80>( unibases, nexts ); break;
  case 88: GetNexts<88>( unibases, nexts ); break;
  case 96: GetNexts<96>( unibases, nexts ); break;
  case 128: GetNexts<128>( unibases, nexts ); break;
  case 144: GetNexts<144>( unibases, nexts ); break;
  case 400: GetNexts<400>( unibases, nexts ); break;
  case 640: GetNexts<640>( unibases, nexts ); break;

  default:
    FatalErr( "GetNexts: Not implemented for K=" << K << "." );
  }
}




#endif
