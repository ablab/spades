///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// SmithWatFree( S, T, gap_penalty )
//
// Let S and T be basevectors or fastavectors, with S very small relative to T.  
// Return the best (lowest) score of an alignment of S with T, relative to the 
// following rules:
//
// (a) a mismatch scores +1
// (b) each blank in a gap scores +2
// (c) each base of S must align with T (so S must be found as a subset of T)
// (d) there is optionally a penalty for gaps on either side of S

// Compile with -funroll-all-loops (perhaps 5% faster).

#include "Alignment.h"
#include "Basevector.h"
#include "CoreTools.h"
#include "Fastavector.h"
#include "ShortVector.h"
#include "dna/Bases.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatFree.h"

template<class L, int MODE>
unsigned int SmithWatFree( const L& S, const L& T, int& best_loc, alignment& a,
			   bool penalize_left_gap, bool penalize_right_gap, 
                           unsigned int mismatch_penalty, unsigned int gap_penalty,
                           unsigned int outer_gap_penalty )

{
     const unsigned int gap_penalty1 = gap_penalty;
     const unsigned int gap_penalty2 = gap_penalty;

     ForceAssertGt( S.size(), 0u );
     ForceAssertGt( T.size(), 0u );

     unsigned int n = S.size( ), N = T.size( );

     ForceAssertLe( n, N );

     avector<char> s, t;
     s.resize(n);
     for ( unsigned int i = 0; i < n; i++ )
          s(i) = S[i];
     t.resize(N);
     for ( unsigned int i = 0; i < N; i++ )
          t(i) = T[i];

     const unsigned int Infinity = 1000000000;
     unsigned int best_score = Infinity;
     avector<unsigned int> x(n+1);
     for ( unsigned int i = 0; i <= n; i++ )
          x(i) = Infinity;
     best_loc = 0;

     for ( unsigned int j = 0; j < N; ++j )
     {    unsigned int* xp = x.x;
          unsigned int lastx = 0;
          if ( penalize_left_gap ) lastx += outer_gap_penalty * j;
	  for ( unsigned int i = 0; i < n; ++i )
	  {    unsigned int ne;
               if ( MODE == 1 ) ne = ( s(i) != t(j) );
               else 
               {    ne = !GeneralizedBase::fromChar( s(i) )
                         .matches( GeneralizedBase::fromChar( t(j) ) );    }
               unsigned int a = lastx + mismatch_penalty * ne;
	       unsigned int b = *xp + gap_penalty2;
	       ++xp;
	       unsigned int c = *xp + gap_penalty1;
	       lastx = *xp;
	       *xp = Min( Min( a, b ), c );    }
	  unsigned int this_score2 = x(n);
	  if ( penalize_right_gap ) this_score2 += outer_gap_penalty * ( N - j - 1 );
          if ( this_score2 <= best_score ) best_loc = j;
          best_score = Min( best_score, this_score2 );    }

     // Redo to find the alignment.

     unsigned int subN = Min( best_loc, 
          (int) (S.size( ) + best_score/(Min(gap_penalty1, gap_penalty2)) + 1) ) + 1;
     t.resize(subN);
     for ( unsigned int i = 0; i < subN; ++i )
          t(i) = T[ best_loc - subN + i + 1 ];
     vec< vec<unsigned char> > from( subN, vec<unsigned char>(n) );
     unsigned int best_score2 = Infinity;
     for ( unsigned int i = 0; i <= n; ++i )
          x(i) = Infinity;
     int best_loc2 = 0;
     for ( unsigned int j = 0; j < subN; ++j )
     {    unsigned int* xp = x.x;
          unsigned int lastx = 0;
	  if ( penalize_left_gap ) lastx += outer_gap_penalty * j;
          for ( unsigned int i = 0; i < n; ++i )
	  {    unsigned int ne;
               if ( MODE == 1 ) ne = ( s(i) != t(j) );
               else 
               {    ne = !GeneralizedBase::fromChar( s(i) )
                         .matches( GeneralizedBase::fromChar( t(j) ) );    }
               unsigned int a = lastx + mismatch_penalty * ne;
               unsigned int b = *xp + gap_penalty2;
               ++xp;
               unsigned int c = *xp + gap_penalty1;
               lastx = *xp;
               if ( a <= b )
               {    if ( a <= c ) from[j][i] = 'a';
                    else from[j][i] = 'c';    }
               else
               {    if ( b <= c ) from[j][i] = 'b';
                    else from[j][i] = 'c';    }
               *xp = Min( Min( a, b ), c );    }
	  unsigned int this_score2 = x(n);
	  if ( penalize_right_gap ) 
          {    this_score2 += outer_gap_penalty 
                    * ((N - best_loc - 1) + (subN - j - 1));    }
          if ( this_score2 <= best_score2 ) best_loc2 = j;
          best_score2 = Min( best_score2, this_score2 );    }
     ForceAssert( best_loc2 == (int) subN-1 );
     int j = best_loc2;
     int i = ((int) n) - 1;
     int lcount = 0, g1count = 0, g2count = 0;
     int last_length = 0;
     avector<int> gaps(0), lengths(0);
     while(1)
     {    if ( from[j][i] == 'a' )
          {    if ( g1count > 0 )
               {    if ( last_length > 0 )
		    {    gaps.Prepend( g1count );
		         lengths.Prepend( last_length );    }
	            g1count = 0;    }
               if ( g2count > 0 )
               {    if ( last_length > 0 )
	            {    gaps.Prepend( -g2count );
                         lengths.Prepend( last_length );    }
                    g2count = 0;    }
               ++lcount;
               --i;
               --j;    }
          else if ( from[j][i] == 'b' )  // gap on long sequence
          {    if ( lcount > 0 )
               {    last_length = lcount;
                    lcount = 0;    }
               ForceAssert( g1count == 0 );
               ++g2count;
               --i;    }
          else                           // gap on short sequence
          {    if ( lcount > 0 )
               {    last_length = lcount;
                    lcount = 0;    }
               ForceAssert( g2count == 0 );
               ++g1count;
               --j;    }
          if ( i < 0 ) break;
          ForceAssert( j >= 0 );    }
     if ( g1count != 0 ) gaps.Prepend( g1count );
     else if ( g2count != 0 ) gaps.Prepend( -g2count );
     else gaps.Prepend(0);
     lengths.Prepend( lcount );
     int pos1 = 0;
     int pos2 = best_loc - subN + 1 + j + 1;
     ForceAssertLe( gaps(0), 0 );
     int first_gap = gaps(0);
     if ( first_gap < 0 )
     {    pos2 -= first_gap;
          gaps(0) = 0;    }
     int errors = best_score;
     a = alignment( pos1, pos2, errors, gaps, lengths );
     ForceAssertEq( best_score, best_score2 );
     return best_score;    }

unsigned int SmithWatFree( const basevector& S, const basevector& T, 
			   int& best_loc, alignment& a,
			   bool penalize_left_gap, bool penalize_right_gap, 
                           unsigned int mismatch_penalty, unsigned int gap_penalty,
                           unsigned int outer_gap_penalty )
{
     return SmithWatFree<basevector,1>( S, T, best_loc, a, penalize_left_gap, 
          penalize_right_gap, mismatch_penalty, gap_penalty, 
          outer_gap_penalty );    }

unsigned int SmithWatFree( const fastavector& S, const fastavector& T, 
			   int& best_loc, alignment& a,
			   bool penalize_left_gap, bool penalize_right_gap, 
                           unsigned int mismatch_penalty, unsigned int gap_penalty,
                           unsigned int outer_gap_penalty )
{
     return SmithWatFree<fastavector,2>( S, T, best_loc, a, penalize_left_gap,
          penalize_right_gap, mismatch_penalty, gap_penalty, 
          outer_gap_penalty );    }
