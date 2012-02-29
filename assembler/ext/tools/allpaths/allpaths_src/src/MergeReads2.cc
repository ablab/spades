///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Alignment.h"
#include "Basevector.h"
#include "CoreTools.h"
#include "math/Functions.h"
#include "MergeReads2.h"
#include "PackAlign.h"
#include "PrintAlignment.h"
#include "Qualvector.h"

// Pick the winning base and compute its score; load the answers into the
// the consensus and quality vectors.  Also keep track of starting positions
// of reads on contigs.

// If a base is called as a gap (so that it doesn't appear in the consensus),
// the scores of the two surrounding bases are lowered by 10.  There must be
// a better approach.

// Note: if scoring = 2, here's how we compute the quality score:
//
//     Computation of scoreA:
//     s1 + ... + sk (sum of scores for agreeing bases)
//     - t1 - ... - td (less scores for disagreeing bases).
//
//     Computation of scoreB:
//     log10(
//          10^s1 + ... + 10^sk (sum of exponentiated scores for agreeing bases)
//          - 10^t1 - ... - 10^td (less exponentiated scores for disagreeing bases)
//          ).   
// *** Note - 4/20/04: what was I thinking??  This might make some sense if we 
// *** computed 10^(s1/10) + ... , and then multiplied by 10 at the end.  As it 
// *** stands this will nearly always return Max(s1,...,sk).  Accordingly, I am 
// *** setting scoreB equal to this Max.
//
//     Computation of reported score: Min( (scoreA+scoreB)/2, 50 ).

// If scoring = 3, using scoring method 2, but then lower the score in accordance
// with the highest conflicting score.

// Modification (10/3/03): scores are first lowered to 50 or less.

void PickWinner

     (  /* inputs */  const vec<unsigned char>& base, 
                      const qualvector& score_orig, const vec<int>& lids,
                      const vec<Bool>& rc, Bool forward, int contig_count,
        /* outputs */ vec<unsigned char>& consensus, qualvector& quality,
                      vec<int>& start_on_contig, vec<int>& stop_on_contig,
                      vec<Bool>& RC, vec<int>& contig_no,
                      int scoring )

{    
     int best_score = 0, best_base = 0;
     const int MinScore = 2, MaxScore = 50;

     // First lower scores so that they are at most 50.

     static qualvector score;
     score.resize( score_orig.size( ) );
     for ( int i = 0; i < (int) score.size( ); i++ )
          score[i] = Min( (unsigned char) 50, score_orig[i] );

     // Set up last_base.

     static unsigned char last_base = NoBase;
     if ( consensus.size( ) == 0 ) last_base = NoBase;

     // Compute scores.

     if ( scoring == 0 || scoring == 1 || scoring == 2 || scoring == 3 )
     {    static vec<int> score_sum(5);
          for ( unsigned int i = 0; i < 5; i++ )
               score_sum[i] = 0;
          for ( unsigned int i = 0; i < base.size( ); i++ )
               score_sum[ base[i] ] += Max( 1, int(score[i]) );
          int max_score = 0;
          for ( unsigned int i = 0; i < 5; i++ )
               if ( score_sum[i] > max_score )
               {    max_score = score_sum[i];
                    best_base = i;    }
          best_score = score_sum[best_base];
          for ( int i = 0; i < 5; i++ )
               if ( i != best_base ) best_score -= score_sum[i];    }

     if ( scoring == 1 || scoring == 2 || scoring == 3)
     {    int m = 0;
          for ( unsigned int i = 0; i < base.size( ); i++ )
               if ( best_base == base[i] ) m = Max( m, int(score[i]) );
          int alt_best_score = m;
          if ( scoring == 1 ) best_score = Min( alt_best_score, best_score );
          else if ( scoring == 2 || scoring == 3 )
          {    int score0 = Min( alt_best_score, best_score );
               int score1 = Max( alt_best_score, best_score );
               best_score = (score0 + score1)/2;    }    }

     if ( scoring == 3 )
     {    int best_conflict = 0;
          for ( int i = 0; i < base.isize( ); i++ )
               if ( base[i] != best_base ) 
                    best_conflict = Max( best_conflict, int(score[i]) );
          best_score = Min( best_score, 2 * (40 - best_conflict) );    }

     // Finalize, all cases.  Takes as input best_base, best_score.

     best_score = Max( best_score, MinScore );
     best_score = Min( best_score, MaxScore );
     if ( last_base == GapBase ) best_score = Max( 0, best_score - 10 );
     if ( best_base != GapBase )
     {    for ( unsigned int i = 0; i < lids.size( ); i++ )
          {    if ( start_on_contig[ lids[i] ] == -1 )
               {    start_on_contig[ lids[i] ] = consensus.size( );
                    RC[ lids[i] ] = (rc[i] == forward);
                    contig_no[ lids[i] ] = contig_count;    }
               stop_on_contig[ lids[i] ] = consensus.size( );    }
          consensus.push_back(best_base);
          quality.push_back(best_score);    }
     else if ( quality.size( ) > 0 ) 
          quality.back( ) = Max( 0, quality.back( ) - 10 );
     last_base = best_base;    }

int Mutations(const basevector& rd1, const basevector& rd2, int p1, int p2, 
     const avector<int>& gaps, const avector<int>& lengths, int nblocks )
{    int answer = 0;
     for ( int j = 0; j < nblocks; j++ )
     {    if ( gaps(j) > 0 ) p2 += gaps(j);
          if ( gaps(j) < 0 ) p1 -= gaps(j);
          for ( int x = 0; x < lengths(j); x++ )
          {    if ( rd1[p1] != rd2[p2] ) ++answer;
               ++p1; ++p2;    }    }
     return answer;    }

// An isolated gap in a pairwise alignment is mobile if it can be shifted one
// base to the left or right, without changing the number of bases which align
// correctly.  The routine CenterMobileGaps is designed to put the gaps into
// "canonical position", so that in a collection of pairwise alignments, the
// gaps which appear in related regions will be in the same position.

// CenterMobileGaps shifts the position of mobile gaps, so that they are as nearly 
// centered as possible, within the range of positions which they can move.
// Furthermore, given a choice between e.g.
//
//                   AGACCTAC         or        AGACCTAC
//                   AGAC-TAC                   AGA-CTAC,
//
// a consistent choice is made.

// Defects in this routine result in serious errors in the output of MergeReads,
// wherein sequences of the form N^n get changed to N^(n+1), where N is any base.

// Cases where an isolated gap could coalesce with another gap are not handled
// correctly.

// Cases where a gap is on the end are not handled correctly.

// The implementation is horrendously inefficient.

Bool CenterMobileGaps( align& a, const basevector& rd1, const basevector& rd2,
     bool verbose, ostream& log )
{    
     int pos1 = a.pos1( ), pos2 = a.pos2( ), nblocks = a.Nblocks( );
     const avector<int> &gaps = a.Gaps( ), &lengths = a.Lengths( );

     int i;
     for ( i = 0; i < nblocks; i++ )
          if ( gaps(i) == 1 || gaps(i) == -1 ) break;
     if ( i == nblocks ) return False;
     Bool changed = False;
     int mutations = Mutations( rd1, rd2, pos1, pos2, gaps, lengths, nblocks );

     avector<int> lengths_orig(0);
     if ( (int) lengths_orig.length < nblocks ) lengths_orig.Setsize(nblocks);
     for ( int j = 0; j < nblocks; j++ )
          lengths_orig(j) = lengths(j);

     int p1 = pos1, p2 = pos2;
     for ( i = 0; i < nblocks; i++ )
     {    if ( gaps(i) == 1 || gaps(i) == -1 )
          {    int left_limit = 0, right_limit = 0;
               if ( i >= 1 )
               {    int p1x = p1, p2x = p2;
                    int save1 = lengths(i-1), save2 = lengths(i);
                    // The middle part of the following for shouldn't be there:
                    for (left_limit = -1; -left_limit <= lengths(i-1); left_limit--)
                    {    a.AddToLength( i-1, -1 );
                         a.AddToLength( i, 1 );
                         if ( gaps(i) == 1 )
                         {    if ( (rd1[p1x-1] == rd2[p2x-1])
                                        != (rd1[p1x-1] == rd2[p2x]) )
                                   break;    }
                         else
                         {    if ( (rd1[p1x-1] == rd2[p2x-1])
                                        != (rd1[p1x] == rd2[p2x-1]) )
                                   break;    }
                         --p1x;
                         --p2x; 
                         if ( -left_limit == lengths(i-1) ) break;    }
                    a.AddToLength( i-1, -left_limit );
                    a.AddToLength( i, left_limit );
                    p1x = p1;
                    p2x = p2;
                    if ( -left_limit != lengths(i-1) ) ++left_limit;
                    // The following two lines are there due to a bug:
                    a.SetLength( i-1, save1 );
                    a.SetLength( i, save2 );
                    // The middle part of the following for shouldn't be there:
                    for ( right_limit = 1; right_limit <= lengths(i); right_limit++ )
                    {    a.AddToLength( i-1, 1 );
                         a.AddToLength( i, -1 );
                         if ( gaps(i) == 1 )
                         {    if ( (rd1[p1x] == rd2[p2x+1])
                                        != (rd1[p1x] == rd2[p2x]) )
                                   break;    }
                         else
                         {    if ( (rd1[p1x+1] == rd2[p2x])
                                        != (rd1[p1x] == rd2[p2x]) )
                                   break;    }
                         ++p1x;
                         ++p2x; 
                         if ( right_limit == lengths(i) ) break;    }
                    a.AddToLength( i-1, -right_limit );
                    a.AddToLength( i, right_limit );
                    if ( right_limit != lengths(i) ) --right_limit;    
                    // The following two lines are there due to a bug:
                    a.SetLength( i-1, save1 );
                    a.SetLength( i, save2 );    }

               int n = (left_limit + right_limit)/2;
               if ( Abs(n) >= 1 )
               {    changed = True;
                    a.AddToLength( i-1, n );
                    a.AddToLength( i, -n );    }

               char the_base = as_base( (gaps(i) == 1) ? rd2[p2] : rd1[p1] );

               if ( left_limit + right_limit - 2*n > 0 
                    && ( the_base == 'A' || the_base == 'C' ) )
               {    changed = True;
                    a.AddToLength( i-1, 1 );
                    a.AddToLength( i, -1 );
                    // The following three lines are there due to a bug:
                    if (Mutations(rd1, rd2, pos1, pos2, gaps, lengths, nblocks) 
                         != mutations)
                    {    a.AddToLength( i-1, -1 );
                         a.AddToLength( i, 1 );    }    }
               if ( left_limit + right_limit - 2*n < 0 
                    && ( the_base == 'G' || the_base == 'T' ) )
               {    changed = True;
                    if ( lengths(i-1) >= 1 ) // shouldn't be needed
                    {    a.AddToLength( i-1, -1 );
                         a.AddToLength( i, 1 );
                         // The following three lines are there due to a bug:
                         if (Mutations(rd1, rd2, pos1, pos2, gaps, lengths, nblocks) 
                              != mutations)
                         {    a.AddToLength( i-1, 1 );
                              a.AddToLength( i, -1 );    }    }    }    }

          if ( gaps(i) > 0 ) p2 += gaps(i);
          if ( gaps(i) < 0 ) p1 -= gaps(i);
          p1 += lengths_orig(i);
          p2 += lengths_orig(i);    }

     if ( changed && !Proper( a, rd1.size( ), rd2.size( ) ) )
     {    if (verbose) log << "compactifying improper alignment between "
               << "reads of lengths " << rd1.size( ) 
               << " and " << rd2.size( ) << endl;
          a.Compactify( rd1.size( ), rd2.size( ) );    }

     if ( changed && verbose ) // XXX
     {    log << "new alignment between reads of lengths " // XXX
               << rd1.size( ) << " and " << rd2.size( ) << endl; // XXX
          PrintVisualAlignment( True, log, rd1, rd2, a );    } // XXX

     return changed;    }
