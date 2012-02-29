///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "ParallelVecUtilities.h"
#include "paths/AlignReadsOnUnibases.h"

void AlignReadsOnUnibases( const vecbvec &jbases,
			   const vecbvec &unibases,
			   vec< triple<int64_t,int,int> > &jaligns,
			   vec<basevector> &jbases_sorted,
			   vec<int64_t> &jbases_sorted_id,
			   ostream &out )
{
     jaligns.clear( );
     jbases_sorted.clear( );
     jbases_sorted_id.clear( );
     int min_read = 1000000000;
     for ( size_t id = 0; id < jbases.size( ); id++ )
          min_read = Min( min_read, jbases[id].isize( ) );
     DPRINT(min_read);
     out << Date( ) << ": sorting jbases" << endl;
     jbases_sorted.resize( jbases.size( ) );
     jbases_sorted_id.resize( jbases.size( ) );
     for ( size_t id = 0; id < jbases.size( ); id++ )
     {    jbases_sorted[id] = jbases[id];
          if ( jbases[id].isize( ) > min_read )
          {    jbases_sorted[id].SetToSubOf( 
                    jbases[id], jbases[id].isize( ) - min_read, min_read );    }
          jbases_sorted_id[id] = id;    }
     ParallelSortSync( jbases_sorted, jbases_sorted_id );
     out << Date( ) << ": looking up unibases" << endl;
     vec< triple<int64_t,int,int> > jaligns0;
     #pragma omp parallel for
     for ( size_t m = 0; m < unibases.size( ); m++ )
     {    vec< triple<int64_t,int,int> > jaligns0m;
          for ( int p = 0; p <= unibases[m].isize( ) - min_read; p++ )
          {    basevector b( unibases[m], p, min_read );
               for ( int pass = 1; pass <= 2; pass++ )
               {    if ( pass == 2 ) b.ReverseComplement( );
                    int64_t low = LowerBound( jbases_sorted, b );
                    int64_t high = UpperBound( jbases_sorted, b );
                    for ( int64_t l = low; l < high; l++ )
                    {    int64_t id = jbases_sorted_id[l];
                         int mx = m;
                         if ( pass == 2 ) mx = -m-1;
                         jaligns0m.push( id, mx, p );    }    }    }
          #pragma omp critical
          {    jaligns0.append(jaligns0m);    }    }
     out << Date( ) << ": sorting aligns" << endl;
     ParallelSort(jaligns0);
     for ( size_t j1 = 0; j1 < jaligns0.size( ); j1++ )
     {    size_t j2;
          for ( j2 = j1 + 1; j2 < jaligns0.size( ); j2++ )
               if ( jaligns0[j2].first != jaligns0[j1].first ) break;
          if ( j2 - j1 == 2 ) 
          {    jaligns.push_back( jaligns0[j1] );
               jaligns.push_back( jaligns0[j1+1] );    }
          j1 = j2 - 1;    }    }

