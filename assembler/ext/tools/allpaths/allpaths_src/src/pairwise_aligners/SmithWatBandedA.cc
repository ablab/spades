///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// #define USE_TIMERS

// SmithWatBandedA( S, T, offset, bandwidth )
//
// Let S and T be basevectors.  Find the best (lowest) scoring alignment of S
// with T, relative to the following rules:
//
// (a) a mismatch scores +1.0;
// (b) each blank in a gap scores +1.5;
// (c) gaps on either end of either read do not count.

// Only search for alignments with the given offset, plus or minus the given
// bandwidth.  The meaning of "offset" is explained by:
//
// example: offset = 2
// S -----------------
// T   ---------------------

// example: offset = -2
// S   ---------------------
// T -----------------

// Return the score of the alignment, and by reference, the alignment itself.

#include "Alignment.h"
#include "Basevector.h"
#include "CoreTools.h"
#include "PackAlign.h"
#include "ShortVector.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatBandedA.h"

template<class X>
float SmithWatBandedA2( const basevector& S, const basevector& T, int offset,
     int bandwidth, align& a, int& errors, ostream *log, int MIS, int INS, int DEL)
{
     const float divider = 2.0;

     // Separately handle the case where the bandwidth is zero, which was not
     // handled correctly in the main code and in any case can be done more
     // efficiently.

     if ( bandwidth == 0 )
     {    errors = 0;
          int start = Max( 0, offset ); 
          int stop = Min( T.isize( ) + offset, S.isize( ) );
          if ( !( start < stop ) ) PRINT3( S.size( ), T.size( ), offset );
          ForceAssertLt( start, stop );
          for ( int i = start; i < stop; i++ )
          {    if ( S[i] != T[ i - offset ] ) errors += MIS;    }
          a.SetNblocks(1);
          a.SetGap( 0, 0 );
          a.SetLength( 0, stop - start );
          a.Setpos1(start);
          a.Setpos2( start - offset );
          return float(errors)/divider;    }

     // Main part of code.
     
     ostream &out = ( log ) ? *log : cout;

     START_TIMER(SW, 4000);
     if ( offset > (int) S.size( ) || offset < - (int) T.size( ) ) {
	 out << "Warning: SmithWatBandedA passed nonsense arguments:" << endl;
	 PRINT2_TO( out, S.size( ), T.size( ) );
	 PRINT2_TO( out, offset, bandwidth );
     }

     // Note that we used to have right = offset + bandwidth.  I added one because
     // there were examples where we didn't get the right alignment.  I'm not sure
     // this is the right fix.

     int left = offset - bandwidth, right = offset + bandwidth + 1;

     const int mismatch_penalty = MIS;
     const int ins_penalty = INS;
     const int del_penalty = DEL;

     const int Infinity = 1000000000;

     int n = (int) S.size( ), N = (int) T.size( );

     vec<char> s;
     s.resize(n);
     for ( int i = 0; i < n; i++ )
          s[i] = S[i];
     int best_score = Infinity;
     int best_j = 0, best_i = 0;
     vec<int> x;
     x.resize(n+1);
     int istart = 0, istop = 0;

     for ( int i = 0; i <= n; i++ )
     {    x[i] = 0;
          if ( !( left <= i && i <= right ) ) x[i] = Infinity;    }

     int jstart = Max( 0, -right-1 ), jstop = Min( N-1, n-left );

     vec< vec<X> > from;
     int from_size = jstop - jstart + 1;
     if ( (int) from.size( ) < from_size ) from.resize( from_size + from_size/5 );
     for ( int i = 0; i < from_size; i++ )
     {    from[i].resize( right - left + 3 );
          unsigned int tempui = from[i].size();
          for ( unsigned int j = 0; j < tempui; j++ )
               from[i][j] = 'n';    }


     // [i,j] to be stored at from[ j - jstart ][ i - left - j + 1 ]

     for ( int j = jstart; j <= jstop; j++ )
     {    x[0] = 0;
          if ( !( left <= -j && -j <= right ) ) x[0] = Infinity;
          int lastx = 0;
          if ( !( left <= -(j-1) && -(j-1) <= right ) ) lastx = Infinity;
          char* sp = &s[0];

          istart = Max( 0, left + j - 1 );
          istop = Min( n - 1, right + j + 1 );
          int* xp = &x[0] + istart;
          if ( istart > 0 )
          {    lastx = *xp;
               *xp = Infinity;    }
          // In SWMIN, got rid of Min ( Min (...) ) since comparisons
          // were already done.

          #define SWMIN(I)                                                         \
               lastx = *xp;                                                        \
               if ( a <= b )                                                       \
               {    if ( a <= c ) {from[j-jstart][I-left-j+1] = 'a'; *xp=a;}       \
                    else {from[j-jstart][I-left-j+1] = 'c'; *xp=c;}    }           \
               else                                                                \
               {    if ( b <= c ) {from[j-jstart][I-left-j+1] = 'b'; *xp=b;}       \
                    else {from[j-jstart][I-left-j+1] = 'c'; *xp=c;}    }

          #define SWCORE(J)                                                        \
               else if ( T[j] == J )                                               \
               {    if ( istart <= istop )                                         \
                    {    int a = lastx + mismatch_penalty * (sp[istart] != J);     \
                         int b = *xp + ins_penalty;                                \
                         ++xp;                                                     \
                         int c = *xp + del_penalty;                                \
                         if ( !( left <= (istart-1) - (j-1)                        \
                              && (istart-1) - (j-1) <= right ) )                   \
                              a = Infinity;                                        \
                         if ( !( istart - (j-1) <= right ) ) b = Infinity;         \
                         if ( !( left <= (istart-1) - j ) ) c = Infinity;          \
                         SWMIN(istart)    }                                        \
                    if ( istart + 1 <= istop )                                     \
                    {    int a = lastx + mismatch_penalty * (sp[istart+1] != J);   \
                         int b = *xp + ins_penalty;                                \
                         ++xp;                                                     \
                         int c = *xp + del_penalty;                                \
                         if ( !( left <= istart - (j-1)                            \
                              && istart - (j-1) <= right ) ) a = Infinity;         \
                         if ( !( istart+1 - (j-1) <= right ) ) b = Infinity;       \
                         if ( !( left <= istart - j ) ) c = Infinity;              \
                         SWMIN(istart+1)    }                                      \
                    for ( int i = istart + 2; i <= istop - 2; i++ )                \
                    {    int a = lastx + mismatch_penalty * (sp[i] != J);          \
                         int b = *xp + ins_penalty;                                \
                         ++xp;                                                     \
                         int c = *xp + del_penalty;                                \
                         SWMIN(i)    }                                             \
                    if ( istop - 1 >= istart + 2 )                                 \
                    {    int a = lastx + mismatch_penalty * (sp[istop-1] != J);    \
                         int b = *xp + ins_penalty;                                \
                         ++xp;                                                     \
                         int c = *xp + del_penalty;                                \
                         if ( !( istop - j <= right ) ) b = Infinity;              \
                         SWMIN(istop-1)    }                                       \
                    if ( istop >= istart + 2 )                                     \
                    {    int a = lastx + mismatch_penalty * (sp[istop] != J);      \
                         int b = *xp + ins_penalty;                                \
                         ++xp;                                                     \
                         int c = *xp + del_penalty;                                \
                         if ( !( istop - j <= right ) ) a = Infinity;              \
                         if ( !( istop - (j-1) <= right ) ) b = Infinity;          \
                         SWMIN(istop)    }    }

          if ( 0 == 1 );
          SWCORE(0)
	  SWCORE(1)
	  SWCORE(2)
	  SWCORE(3)

          if ( istop < n - 1 )
          {    ++xp;
               *xp = Infinity;    }

          if ( istop == n - 1 ) 
          {    if ( x[n] < best_score ) 
               {    best_j = j;
                    best_i = n-1;    }
               best_score = Min( best_score, x[n] );    }    }

     for ( int i = istart; i <= istop; i++ )
     {    if ( x[i] < best_score )
          {    
               // cout << "type 2 set\n"; // XXX
               // PRINT3( i, istart, istop ); // XXX
               // PRINT( x[i] ); // XXX
               best_j = jstop;
               best_i = i - 1;     // NEW NEW NEW NEW NEW NEW!!!!!
               // best_i = i;    
               // PRINT2( best_i, best_j ); // XXX
                    }
          best_score = Min( x[i], best_score );    }

     ForceAssert( best_i < (int) S.size( ) );
     ForceAssert( best_j < (int) T.size( ) );

     int j = best_j, i = best_i;
     int lcount = 0, g1count = 0, g2count = 0, last_length = 0;
     vec<int> gaps, lengths;
     gaps.resize(0);
     lengths.resize(0);

     // int a_count = 0; // XXX
     while(1)
     {    
          if ( j < 0 || i < 0 || j-jstart < 0 || i-left-j+1 < 0 
               || i-left-j+1 >= right-left+3 ) 
               break;

          /*
          char letter = from[j-jstart][i-left-j+1]; // XXX
          if ( letter == 'a' ) ++a_count; // XXX
          else  // XXX
          {    if ( a_count > 0 ) cout << " a^" << a_count << " "; // XXX
               a_count = 0; // XXX
               cout << letter;    } // XXX
          */

          if ( from[j-jstart][i-left-j+1] == 'a' )
          {    if ( g1count > 0 )
               {    gaps.push_back( g1count );
                    lengths.push_back( last_length );
                    g1count = 0;    }
               if ( g2count > 0 )
               {    gaps.push_back( -g2count );
                    lengths.push_back( last_length );
                    g2count = 0;    }
               if ( i == 0 || j == 0 ) break; // NEW NEW NEW !!!
               ++lcount;
               --i;
               --j;    }
          else if ( from[j-jstart][i-left-j+1] == 'b' )  // gap on long sequence
          {    if ( lcount > 0 )
               {    last_length = lcount;
                    lcount = 0;    }
               ForceAssert( g1count == 0 );
               ++g2count;
               --i;    }
          else if ( from[j-jstart][i-left-j+1] == 'c' )  // gap on short sequence
          {    if ( lcount > 0 )
               {    last_length = lcount;
                    lcount = 0;    }
               ForceAssert( g2count == 0 );
               ++g1count;
               if ( j == 0 ) break; // NEW NEW NEW !!!
               --j;    }
          else break;    }
          // if ( i < 1 ) break;
          // ForceAssert( j >= 0 );    }
     // if ( a_count > 0 ) cout << " a^" << a_count << " "; // XXX
     // cout << "\n"; // XXX
 
     // ForceAssert( g1count == 0 );
     // ForceAssert( g2count == 0 );
     // gaps.push_back(0);

     if ( g1count != 0 ) gaps.push_back( g1count );
     else if ( g2count != 0 ) gaps.push_back( -g2count );
     else gaps.push_back(0);

     lengths.push_back( lcount + 1 );
     // lengths.push_back(lcount);

     int pos1 = i, pos2 = j;
     if ( pos1 < 0 || pos2 < 0 )
     {    ++pos1;
          ++pos2;    }

     errors = best_score;
     a.Setpos1(pos1);
     a.Setpos2(pos2);
     a.SetNblocks( gaps.size( ) );
     int nb = gaps.size( );
     for ( int i = 0; i < nb; i++ )
     {    a.SetGap( nb - i - 1, gaps[i] );
          a.SetLength( nb - i - 1, lengths[i] );    }

     if ( pos1 < 0 || pos2 < 0 || a.Pos1( ) > (int) S.size( )
          || a.Pos2( ) > (int) T.size( ) )
     {    
          // See no evil...

          /*
          cout << "Warning: SmithWatBandedA produced a garbage alignment, "
               << "presumably because of a bug.\n";
          cout << "Replacing the garbage alignment by an alignment which has "
               << "one-base overlap.\n";
          PRINT3( pos1, a.Pos1( ), S.size( ) );
          PRINT3( pos2, a.Pos2( ), T.size( ) );
          int pos1x, pos2x, errors;
          avector<int> gapsx, lengthsx;
          a.Unpack( pos1x, pos2x, errors, gapsx, lengthsx );
          cout << "gaps/lengths:";
          for ( unsigned int i = 0; i < lengthsx.length; i++ )
               cout << " " << gapsx(i) << "/" << lengthsx(i);
          cout << endl;
          PRINT2(offset, bandwidth);
          if ( S.size( ) < 1000 ) PRINT( S.ToString( ) );
          if ( T.size( ) < 1000 ) PRINT( T.ToString( ) );
          */

          ForceAssert( S.size( ) > 0 );
          ForceAssert( T.size( ) > 0 );
          Bool base_matches = ( S[0] == T[ T.size( ) - 1 ] );

          avector<int> gapsz(1), lengthsz(1);
          gapsz(0) = 0;
          lengthsz(0) = 1;
          errors = (base_matches ? 0 : 1);
          a.Set( 0, (int) T.size( ) - 1, gapsz, lengthsz );

          STOP_TIMER(SW);
          if (base_matches) return 0;
          else return float(mismatch_penalty)/divider;    }

     STOP_TIMER(SW);
     return float(best_score)/divider;    }

template float SmithWatBandedA2<unsigned char>( const basevector&, 
     const basevector&, int, int, align&, int&, ostream*, int, int, int );
template float SmithWatBandedA2<unsigned short>( const basevector&, 
     const basevector&, int, int, align&, int&, ostream*, int, int, int );
template float SmithWatBandedA2<unsigned int>( const basevector&, 
     const basevector&, int, int, align&, int&, ostream*, int, int, int );
