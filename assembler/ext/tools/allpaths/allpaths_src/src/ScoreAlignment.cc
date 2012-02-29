///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <math.h>

#include "Alignment.h"
#include "math/Arith.h"
#include "Basevector.h"
#include "math/Functions.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "ScoreAlignment.h"

// These functors exist to prevent duplication of code to handle rc'ed
// elements.  If the qual() method is handed an empty qualvector, it
// presumes the base is of infinite (255) quality.

struct null_coord_xform 
{
  char base( const basevector& bv, const int coord ) const
    { return bv[ coord ]; }
  unsigned char qual( const qualvector& qv, const int coord ) const
    { return ( qv.size() == 0 ? 255 : qv[ coord ] ); }
};
  
struct rc_coord_xform 
{
  char base( const basevector& bv, const int coord ) const
    { return 3 - bv[ bv.size() - coord - 1]; }
  unsigned char qual( const qualvector& qv, const int coord ) const
    { return ( qv.size() == 0 ? 255 : qv[ qv.size() - coord - 1 ] ); }
};
  

// Note.  The routines in this file penalize errors in alignments 
// "non-logarithmically".  Thus, a single error at a score of quality 40 is as
// bad as 10 errors at scores of quality 30.

// A score of 255 is considered perfect.

Float Prob( unsigned char q )
{    static Bool first_pass(True);
     static Float prob[256];
     if (first_pass)
     {    for ( int i = 0; i < 256; i++ )
               prob[i] = Pow( 10.0, Float(i) / Float(-10.0) );
          prob[255] = 0.0;
          first_pass = False;    }
     return prob[q];    }

// A match has a score of 0.  A mismatch is scored in the following
// way: For each sequence, we find the error rate (i.e., the
// probability of being in error) of the base with the lowest quality
// score among the mismatched base and the bases on either side of it.
// The score for the match is the reciprocal of the sum of these two
// probabilities.

template< class rd2_transform >
Float ScoreMatch( unsigned int p1, unsigned int p2,
		  const basevector& rd1, const qualvector& scores1, 
		  const basevector& rd2, const qualvector& scores2,
		  const rd2_transform& trans )
{
  if ( rd1[p1] == trans.base(rd2,p2) ) 
    return 0;

  Float pr1, pr2;
  pr1 = Prob(scores1[p1]);
  if ( p1 > 0 )
    pr1 = Max( pr1, Prob(scores1[p1-1]) );
  if ( p1 < (scores1.size() - 1) )
    pr1 = Max( pr1, Prob(scores1[p1+1]) );
	       
  pr2 = Prob(trans.qual(scores2,p2));
  if ( p2 > 0 )
    pr2 = Max( pr2, Prob(trans.qual(scores2,p2-1) ) );
  if ( p2 < (scores2.size() - 1) )
    pr2 = Max( pr2, Prob(trans.qual(scores2,p2+1) ) );
  
  return Float(1) / (pr1+pr2);    
}

// A gap is scored as follows: First, the probabilities that the bases
// on either side of the gap (the flanking bases) are summed. Then we
// find the lowest quality base in the other sequence that corresponds
// to either the gap or one of the flanking bases.  The probability of
// this base is added to the sum, we take the reciprocal and multiply
// that by twice the length of the gap.

// Note that this fuzziness means that a mobile gap may not end up in
// the best place:
//
// 4044                             4044
// 0 00                             0 00
//
// ATTG                             ATTG
//  |    has a score of ~2, as does   |
// A TG                             AT G
//
// 4 44                             44 4
// 0 00                             00 0
//
// even though the first alignment could be considered better.

template< class rd1_transform, class rd2_transform >
Float ScoreGapOnFirstRead( int p1, int p2, int gap,
			   const basevector& rd1, const qualvector& scores1, 
			   const basevector& rd2, const qualvector& scores2,
			   const rd1_transform& rd1_trans, 
			   const rd2_transform& rd2_trans )
{
  if ( gap == 0 ) 
    return 0;
  
  Float pr1 = Prob( rd1_trans.qual( scores1, p1 ) );

  if ( p1 > 0 )
    pr1 += Prob( rd1_trans.qual( scores1, p1-1 ) );

  Float pr2 = 0;

  int gap_begin = Max( p2 - 1, 0 );
  int gap_end = Min( p2 + gap + 1, (int) rd2.size() );
  for ( int gap_base = gap_begin; gap_base < gap_end; ++gap_base )
    pr2 = Max( pr2, Prob( rd2_trans.qual( scores2, gap_base ) ) );
  
  // If both sequences are "perfect" (i.e. they have no quality
  // scores), both probabilities will be zero, which will cause an
  // arithmetic exception in the return statement.  To avoid this, we
  // set the sum of their probabilities to that of a disagreement
  // between two bases of actual quality 255, i.e. something
  // fantastically small.

  Float sum = pr1+pr2;

  if ( sum == 0.0 )
      sum = Pow( 10.0, (255.0 + 255.0) / -10.0 );

  return Float(2) * Float(gap) * Float(1) / ( sum );
}

template< class rd1_transform, class rd2_transform >
Float ScoreGapOnSecondRead( int p1, int p2, int gap,
			    const basevector& rd1, const qualvector& scores1, 
			    const basevector& rd2, const qualvector& scores2,
			    const rd1_transform& rd1_trans, 
			    const rd2_transform& rd2_trans )
{
  return ScoreGapOnFirstRead( p2, p1, gap, 
			      rd2, scores2,
			      rd1, scores1,
			      rd2_trans,
			      rd1_trans );
}

// Given an alignment between two sequences, plus quality information for both
// sequences or for the first sequence, score the alignment.

template< class rd2_transform >
Float ScoreAlignment( const align& a, 
		      const basevector& rd1, const qualvector& scores1, 
		      const basevector& rd2, const qualvector& scores2,
		      const rd2_transform& trans, int start1, int stop1,
                      int start2, int stop2, Bool ignore_gaps )

{    
     int pos1 = a.pos1( ), pos2 = a.pos2( );
     if ( stop1 < 0 ) stop1 = rd1.size( );
     if ( stop2 < 0 ) stop2 = rd2.size( );
     ForceAssertGe( pos1, 0 );
     ForceAssertGe( pos2, 0 );
     const avector<int>& gaps = a.Gaps( );
     const avector<int>& lengths = a.Lengths( );
     int nblocks = a.Nblocks( );

     Float score = 0; // bigger is worse

     int p1 = pos1, p2 = pos2, total_len = 0;

     for ( int j = 0; j < nblocks; j++ )
     {    
          // We ignore gaps at the beginning and end.

          // ***** THIS IS DEFECTIVE CODE.  The lines p2 += gaps(j) and      *****
          // ***** p1 -= gaps(j) should always be run -- whether the gap is  *****
          // ***** internal or not.  However, there may never be gaps at the *****
          // ***** beginning or end.  REGARDLESS, THIS CODE SHOULD BE FIXED. *****

          if ( gaps(j) != 0 && j > 0 && (j < nblocks-1 || lengths(j) > 0) )
          {
               if ( gaps(j) > 0 ) 
               {    if ( !ignore_gaps && p1 >= start1 && p1 < stop1 
                         && p2 >= start2 && p2 < stop2 )
                         score += ScoreGapOnFirstRead( p1, p2, gaps(j),
						       rd1, scores1, rd2, scores2, 
						       null_coord_xform(), trans );
                    p2 += gaps(j);    }

               else if ( gaps(j) < 0 )
               {    if ( !ignore_gaps && p1 >= start1 && p1 < stop1 
                         && p2 >= start2 && p2 < stop2 )
                         score += ScoreGapOnSecondRead( p1, p2, -gaps(j),
						        rd1, scores1, rd2, scores2, 
						        null_coord_xform(), trans );
                    p1 -= gaps(j);    }    }

          for ( int i = 0; i < lengths(j); i++ )
          {    if (p1+i >= start1 && p1+i < stop1 && p2+i >= start2 && p2+i < stop2)
               {      Float match_score 
                         = ScoreMatch( p1+i, p2+i, rd1, scores1, rd2, scores2,
				       trans );
                      score += match_score;    }    }

          p1 += lengths(j);
          p2 += lengths(j);

          total_len += lengths(j);    }

     // The following weird maneuvering overcomes some floating point
     // reproducibility problems.  I'm not sure how many of the intermediate
     // steps are actually needed.

     static Float a1, a1b, a2, answer;
     a1 = Pow( Float(total_len) + Float(0.0000001), Float(1.5) );
     a1b = Float(1) + score;
     a2 = Float(1000) * a1b;
     answer = a2/a1;
     return answer;    }

Float ScoreAlignment( const align& a, 
		      const basevector& rd1, const qualvector& scores1, 
		      const basevector& rd2, const qualvector& scores2,
                      int start1, int stop1, int start2, int stop2,
                      Bool ignore_gaps )
{
  return ScoreAlignment( a, rd1, scores1, rd2, scores2, null_coord_xform(),
                         start1, stop1, start2, stop2, ignore_gaps );
}
			 
Float ScoreAlignment( Bool rd2_is_rc, const align& a, 
		      const basevector& rd1, const qualvector& scores1, 
		      const basevector& rd2, const qualvector& scores2,
                      int start1, int stop1, int start2, int stop2,
                      Bool ignore_gaps )
{    
  if (rd2_is_rc) 
    return ScoreAlignment( a, rd1, scores1, rd2, scores2, rc_coord_xform(),
                           start1, stop1, start2, stop2, ignore_gaps );
  return ScoreAlignment( a, rd1, scores1, rd2, scores2, null_coord_xform(),
                         start1, stop1, start2, stop2, ignore_gaps );   
}

// ScoreAlignmentPoly: Like ScoreAlignment, but different in the following ways:
// 1. Gaps are not penalized.
// 2. Instead of returning a score like ScoreAlignment would, we return the number
//    of mismatches that would have to be ignored in order to lower the 
//    ScoreAlignment-type score to <= 50.
// 3. We substitute (sum of quality scores)/40 for total_length.

template< class rd2_transform >
int ScoreAlignmentPoly( const align& a, 
		      const basevector& rd1, const qualvector& scores1, 
		      const basevector& rd2, const qualvector& scores2,
		      const rd2_transform& trans, int start1, int stop1,
                      int start2, int stop2 )

{    int pos1 = a.pos1( ), pos2 = a.pos2( );
     if ( stop1 < 0 ) stop1 = rd1.size( );
     if ( stop2 < 0 ) stop2 = rd2.size( );
     ForceAssertGe( pos1, 0 );
     ForceAssertGe( pos2, 0 );
     const avector<int> &gaps = a.Gaps( ), &lengths = a.Lengths( );
     int nblocks = a.Nblocks( );

     int p1 = pos1, p2 = pos2, total_qual = 0;

     static vec<Float> mismatches;
     mismatches.clear( );
     for ( int j = 0; j < nblocks; j++ )
     {    
          if ( gaps(j) > 0 ) p2 += gaps(j);
          else if ( gaps(j) < 0 ) p1 -= gaps(j);

          for ( int i = 0; i < lengths(j); i++ )
          {    if (p1+i >= start1 && p1+i < stop1 && p2+i >= start2 && p2+i < stop2)
               {    Float match_score = 
                         ScoreMatch( p1+i, p2+i, rd1, scores1, rd2, scores2, trans );
                    if ( match_score > Float(0.0) ) mismatches.push_back(match_score);    }
               if ( scores2.size( ) == 0 ) total_qual += scores1[p1+i];
               else total_qual += ( scores1[p1+i] + scores2[p2+i] ) / 2;    }

          p1 += lengths(j);
          p2 += lengths(j);    }

     sort( mismatches.begin( ), mismatches.end( ) );

     Float length_pow 
          = Pow( Float(total_qual)/Float(40) + Float(0.0000001), Float(1.5) );
     Float score = 0.0;
     int i;
     for ( i = 0; i < (int) mismatches.size( ); i++ )
     {    score += mismatches[i];
          if ( Float(1000) * Float(Float(1)+score)/length_pow > Float(50) ) 
               break;    }
     return (int) mismatches.size( ) - i;    }

int ScoreAlignmentPoly( const align& a, 
		      const basevector& rd1, const qualvector& scores1, 
		      const basevector& rd2, const qualvector& scores2,
                      int start1, int stop1, int start2, int stop2 )
{
  return ScoreAlignmentPoly( a, rd1, scores1, rd2, scores2, null_coord_xform(),
      start1, stop1, start2, stop2 );
}
			 
int ScoreAlignmentPoly( Bool rd2_is_rc, const align& a, 
		      const basevector& rd1, const qualvector& scores1, 
		      const basevector& rd2, const qualvector& scores2,
                      int start1, int stop1, int start2, int stop2 )
{    
  if (rd2_is_rc) 
    return ScoreAlignmentPoly( a, rd1, scores1, rd2, scores2, rc_coord_xform(),
         start1, stop1, start2, stop2 );
  return ScoreAlignmentPoly( a, rd1, scores1, rd2, scores2, null_coord_xform(),
         start1, stop1, start2, stop2 );
}

// Regap shifts the position of gaps in an alignment so as to improve the
// quality of the alignment, relative to given quality scores.

// At present, only isolated (single-space) gaps are handled.  The procedure is
// as follows.  We say that a gap is left-mobile if when "exchanged" with the
// base on its left, the naive score of the alignment (computed without reference 
// to quality scores) does not worsen.  Similarly we can define a right-mobile
// gap.  Given this framework, what we do is find (amongst all the positions
// that a gap can move to) that position which gives the best score, relative to
// the quality scores.

template< class rd2_transform >
void Regap( align& a, 
	    const basevector& rd1, const qualvector& scores1,
	    const basevector& rd2, const qualvector& scores2,
	    const rd2_transform& trans )

{    Bool changed = False;
     Bool zero_block_created = False;

     int pos1 = a.pos1( ), pos2 = a.pos2( );
     int nblocks = a.Nblocks( );

     int p1 = pos1, p2 = pos2;
     
     const avector<int>& gaps = a.Gaps();
     const avector<int>& lengths = a.Lengths();

     for ( int j = 0; j < nblocks; j++ )
     {    
       // We ignore gaps at the beginning and end.
       
       if ( gaps(j) != 0 && j > 0 && (j < nblocks-1 || lengths(j) > 0) )
       {
	 if ( gaps(j) > 0 ) 
	 {    // There's a gap of gaps(j) on read 1.
	   if ( gaps(j) == 1 )
	   {    // Try to find a better place to the left for the gap.
	     int shift = 0, best_shift = 0;
	     Float score = ScoreGapOnFirstRead( p1, p2, 1, 
						rd1, scores1, rd2, scores2, 
						null_coord_xform(), trans );
	     Float best_score = score;
	     while( --shift >= -lengths(j-1) )
	     {    
	       int o = rd1[p1+shift] == trans.base( rd2, p2+shift );
	       int n = rd1[p1+shift] == trans.base( rd2, p2+shift+1 );
	       if ( o && !n ) break;
	       score -= ScoreGapOnFirstRead( p1+shift+1, p2+shift+1, 1, 
					     rd1, scores1, rd2, scores2,
					     null_coord_xform(), trans );
	       score += ScoreGapOnFirstRead( p1+shift, p2+shift, 1, 
					     rd1, scores1, rd2, scores2,
					     null_coord_xform(), trans );
	       if ( !o )
		 score -= ScoreMatch( p1+shift, p2+shift,
				      rd1, scores1, rd2, scores2,
				      trans );
	       if ( !n )
		 score += ScoreMatch( p1+shift, p2+shift+1,
				      rd1, scores1, rd2, scores2,
				      trans );
	       if ( score < best_score - Float(0.01) )
	       {
		 best_score = score;
		 best_shift = shift;   
	       }
	     }
	     // Try to find a better place to the right for the gap.
	     score = ScoreGapOnFirstRead( p1, p2, 1, 
					  rd1, scores1, rd2, scores2, 
					  null_coord_xform(), trans );
	     shift = 0;
	     while( ++shift <= lengths(j) )
	     {
	       if ( p1+shift >= (int) rd1.size( ) ||
		    p2+shift >= (int) rd2.size( ) ) break;
	       int o = rd1[p1+shift-1] == trans.base( rd2, p2+shift );
	       int n = rd1[p1+shift-1] == trans.base( rd2, p2+shift-1 );
	       if ( o && !n ) break;
	       score -= ScoreGapOnFirstRead( p1+shift-1, p2+shift-1, 1, 
					     rd1, scores1, rd2, scores2,
					     null_coord_xform(), trans );
	       score += ScoreGapOnFirstRead( p1+shift, p2+shift, 1, 
					     rd1, scores1, rd2, scores2,
					     null_coord_xform(), trans );
	       if ( !o )
		 score -= ScoreMatch( p1+shift-1, p2+shift,
				      rd1, scores1, rd2, scores2, 
				      trans );
	       if ( !n )
		 score += ScoreMatch( p1+shift-1, p2+shift-1,
				      rd1, scores1, rd2, scores2,
				      trans );
	       if ( score < best_score - Float(0.01) )
	       {
		 best_score = score;
		 best_shift = shift;
	       }
	     }
	     if ( best_shift != 0 )
	     {
	       // Shift the gap.
	       changed = True;
	       a.AddToLength( j-1, best_shift );
	       a.AddToLength( j, -best_shift );
	       p1 += best_shift;
	       p2 += best_shift;
	       if ( lengths(j-1) == 0 || lengths(j) == 0 )
		 zero_block_created = True;
	     }
	   }
	   p2 += gaps(j);    
	 }
	 
	 else if ( gaps(j) < 0 )
	 {    // There's a gap of -gaps(j) on read 2.
	   if ( gaps(j) == -1 )
	   {    // Try to find a better place to the left for the gap.
	     int shift = 0, best_shift = 0;
	     Float score = ScoreGapOnSecondRead( p1, p2, 1, 
						 rd1, scores1, rd2, scores2,
						 null_coord_xform(), trans );
	     Float best_score = score;
	     while( --shift >= -lengths(j-1) )
	     {    
	       int o = rd2[p2+shift] == trans.base( rd1, p1+shift );
	       int n = rd2[p2+shift] == trans.base( rd1, p1+shift+1 );
	       if ( o && !n ) break;
	       score -= ScoreGapOnSecondRead( p1+shift+1, p2+shift+1, 1, 
					      rd1, scores1, rd2, scores2,
					      null_coord_xform(), trans );
	       score += ScoreGapOnSecondRead( p1+shift, p2+shift, 1, 
					      rd1, scores1, rd2, scores2,
					      null_coord_xform(), trans );
	       if ( !o )
		 score -= ScoreMatch( p1+shift, p2+shift,
				      rd1, scores1, rd2, scores2,
				      trans );
	       if ( !n )
		 score += ScoreMatch( p1+shift+1, p2+shift,
				      rd1, scores1, rd2, scores2,
				      trans );
	       if ( score < best_score - Float(0.01) )
	       {
		 best_score = score;
		 best_shift = shift;  
	       }
	     }
	     // Try to find a better place to the right for the gap.
	     score = ScoreGapOnSecondRead( p1, p2, 1, 
					   rd1, scores1, rd2, scores2, 
					   null_coord_xform(), trans );
	     shift = 0;
	     while( ++shift <= lengths(j) )
	     {
	       if ( p1+shift >= (int) rd1.size( ) ||
		    p2+shift >= (int) rd2.size( ) ) break;
	       int o = rd2[p2+shift-1] == trans.base( rd1, p1+shift );
	       int n = rd2[p2+shift-1] == trans.base( rd1, p1+shift-1 );
	       if ( o && !n ) break;
	       score -= ScoreGapOnSecondRead( p1+shift-1, p2+shift-1, 1, 
					      rd1, scores1, rd2, scores2,
					      null_coord_xform(), trans );
	       score += ScoreGapOnSecondRead( p1+shift, p2+shift, 1, 
					      rd1, scores1, rd2, scores2,
					      null_coord_xform(), trans );
	       if ( !o )
		 score -= ScoreMatch( p1+shift, p2+shift-1,
				      rd1, scores1, rd2, scores2,
				      trans );
	       if ( !n )
		 score += ScoreMatch( p1+shift-1, p2+shift-1,
				      rd1, scores1, rd2, scores2,
				      trans );
	       if ( score < best_score - Float(0.01) )
	       {
		 best_score = score;
		 best_shift = shift;    
	       }
	     }
	     if ( best_shift != 0 )
	     {    
	       // Shift the gap.
	       changed = True;
	       a.AddToLength( j-1, best_shift );
	       a.AddToLength( j, -best_shift );
	       p1 += best_shift;
	       p2 += best_shift;
	       if ( lengths(j-1) == 0 || lengths(j) == 0 )
		 zero_block_created = True;   
	     }
	   }
	   p1 -= gaps(j);  
	 }
       }
       
       p1 += lengths(j);
       p2 += lengths(j);    
     }
     
     if ( changed && zero_block_created )
       a.Compactify( rd1.size( ), rd2.size( ) );    
} 


void Regap( align& a, 
	    const basevector& rd1, const qualvector& scores1, 
	    const basevector& rd2, const qualvector& scores2 )
{
  Regap( a, rd1, scores1, rd2, scores2, null_coord_xform() );    
}

void Regap( Bool rd2_is_rc, 
	    align& a, 
	    const basevector& rd1, const qualvector& scores1, 
	    const basevector& rd2, const qualvector& scores2 )
{
  if (rd2_is_rc) 
    Regap( a, rd1, scores1, rd2, scores2, rc_coord_xform() );
  else 
    Regap( a, rd1, scores1, rd2, scores2, null_coord_xform() );    
}
