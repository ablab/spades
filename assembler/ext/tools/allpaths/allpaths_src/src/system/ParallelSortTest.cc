///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ParallelSortTest.  Time the sorting of a large vector.
//
// Some results for our hardware:
//
// arguments       date      machine    memory to use    time used in sort
//                                      (GB)             (minutes)
//
// X=4  D=1 P=32   5/18/11   crd4        32              0.61      
// X=16 D=1 P=32   5/18/11   crd4       128              2.59      
//
// (Note that running the program will take substantially longer than the sort.)

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
#define _GLIBCXX_PARALLEL

#include <omp.h>

#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "system/SortInPlace.h"

template<int D> class blob {

     public:

     int x[D];

     void SetFunny( int m )
     {    m += 384794311;
          for ( int i = 0; i < D; i++ )
          {    x[i] = (m*m) % 123456789;
               if ( i > 0 ) x[i] += ( x[i-1] * x[i-1] ) % 987654321;    }    }

     friend Bool operator<( const blob& b1, const blob& b2 )
     {    for ( int j = 0; j < D; j++ )
          {    if ( b1.x[j] < b2.x[j] ) return True;
               if ( b1.x[j] > b2.x[j] ) return False;    }
          return False;    }

     friend int compare( blob const& b1, blob const& b2 )
     {
         int result = 0;
         int const* itr = b1.x;
         int const* end = b1.x+D;
         int const* itr2 = b2.x;
         for ( ; !result && itr != end; ++itr, ++itr2 )
             result = compare(*itr,*itr2);
         return result;
     }
};

static Bool In_place_sort, Parallel_sort;

template<int D> void Slobber( const longlong N )
{
    vec<blob<D> > x(N);
    for ( longlong i = 0; i < N; i++ )
        x[i].SetFunny(i);
    vec<blob<D> > y(x);

    double clock;

    if (In_place_sort)
    {    clock = WallClockTime();
         sortInPlaceParallel(x.begin(),x.end());
         cout << TimeSince(clock) << " used in InPlaceSort" << endl;
         for ( longlong i = 1; i < N; i++ )
         {    if ( x[i] < x[i - 1] )
              {    cout << "Sort failed" << endl;
                   break;    }    }    }

    if (Parallel_sort)
    {    clock = WallClockTime();
         ParallelSort(y);
         cout << TimeSince(clock) << " used in ParallelSort" << endl;
         for ( longlong i = 1; i < N; i++ )
         {    if ( y[i] < y[i - 1] )
              {    cout << "ParallelSort failed" << endl;
                   break;    }    }    }    }

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Double_Doc(X, "number of entries in vector, in billions");
     CommandArgument_Int_OrDefault_Doc(D, 1, "size of data structure in ints");
     CommandArgument_Int_OrDefault_Doc(P, 0, "number of processors");
     CommandArgument_Bool_OrDefault(IN_PLACE_SORT, False);
     CommandArgument_Bool_OrDefault(PARALLEL_SORT, True);
     EndCommandArguments;

     if ( P == 0 ) P = omp_get_max_threads( );

     In_place_sort = IN_PLACE_SORT;
     Parallel_sort = PARALLEL_SORT;

     cout << "memory to use = " 
          << X * D * sizeof(int) * ( PARALLEL_SORT ? 2 : 1 ) << " GB" << endl;

     omp_set_num_threads(P);
     longlong N = longlong( round( X * 1000.0 * 1000.0 * 1000.0 ) );
     if      ( D == 1 )  Slobber<1>(N);
     else if ( D == 2 )  Slobber<2>(N);
     else if ( D == 4 )  Slobber<4>(N);
     else if ( D == 6 )  Slobber<6>(N);
     else if ( D == 8 )  Slobber<8>(N);
     else if ( D == 10 ) Slobber<10>(N);
     else if ( D == 12 ) Slobber<12>(N);
     else if ( D == 14 ) Slobber<14>(N);
     else if ( D == 16 ) Slobber<16>(N);
     else
     {    cout << "Not implemented for that value of D." << endl;
          exit(1);    }    }
