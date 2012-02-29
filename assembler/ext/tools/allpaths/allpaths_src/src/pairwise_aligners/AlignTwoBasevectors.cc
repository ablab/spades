// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// #define USE_TIMERS

// MergeTwoBasevectors( b1, b2, a, c, q1, q2, q ):
//
// Given two basevectors b1, b2 and an alignment a of them, produce a merged
// basevector c.  Also merge quality scores q1, q2 to yield q.

// This is dog-food code, but under optimal conditions, it will do no 
// irreversible damage to the world.

// The code here is adapted from code in PrintVisualAlignment.

#include "Alignment.h"
#include "pairwise_aligners/AlignTwoBasevectors.h"
#include "math/Arith.h"
#include "Badness.h"
#include "Basevector.h"
#include "kmers/KmerShape.h"
#include "CoreTools.h"
#include "math/Functions.h"
#include "pairwise_aligners/MakeAlignsStatic.h"
#include "Overlap.h"
#include "PackAlign.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "ScoreAlignment.h"
#include "ShortVector.h"

// #include <time.h> // AAA

void MergeTwoBaseVectors( const basevector& b1, const basevector& b2,
     const align& a, basevector& c, const qualvector& q1,
     const qualvector& q2, qualvector& q )
{    int pos1 = a.pos1( ), pos2 = a.pos2( ), Pos1 = a.Pos1( ), Pos2 = a.Pos2( );
     const avector<int> &gaps = a.Gaps( ), &lengths = a.Lengths( );
     ForceAssertEq( b1.size( ), q1.size( ) );
     ForceAssertEq( b2.size( ), q2.size( ) );
     ForceAssertLt( pos2, (int) b2.size( ) );
     int nblocks = a.Nblocks( );

     int errors = ActualErrors( b1, b2, a );

     c.Setsize( b1.size( ) + b2.size( ) + errors );
     q.resize( b1.size( ) + b2.size( ) + errors );

     int i = 0, p1 = pos1, p2 = pos2;

     if ( pos1 > 0 && pos2 > 0 )
     {    cout << "MergeTwoBasevectors: I've got pos1 = " << pos1 <<
               " and pos2 = " << pos2 << ", which is not what I expected.\n";
          ForceAssert( 0 == 1 );    }

     if ( pos2 == 0 )
          for ( ; i < pos1; i++ )
          {    c.Set( i, b1[i] );
               q[i] = q1[i];    }
     else for ( ; i < pos2; i++ )
          {    c.Set( i, b2[i] );
               q[i] = q2[i];    }
     for ( int j = 0; j < nblocks; j++ )
     {    if ( gaps(j) > 0 )
               for ( int x = 0; x < gaps(j); x++ )
               {    c.Set( i, b2[p2] );
                    q[i] = q2[p2];
                    //q[i] = 0;
                    ++p2;
                    ++i;    }
          if ( gaps(j) < 0 )
               for ( int x = 0; x < -gaps(j); x++ )
               {    c.Set( i, b1[p1] );
	            q[i] = q1[p1];
                    //q[i] = 0;
                    ++p1;
                    ++i;    }
          for ( int x = 0; x < lengths(j); x++ )
          {    if ( b1[p1] == b2[p2] )
               {    c.Set( i, b1[p1] );
                    q[i] = Max( q1[p1], q2[p2] );    }
               else if ( q1[p1] >= q2[p2] )
               {    c.Set( i, b1[p1] );
                    q[i] = q1[p1] - q2[p2];    }
               else
               {    c.Set( i, b2[p2] );
                    q[i] = q2[p2] - q1[p1];    }
               ++p1;
               ++p2;
               ++i;    }    }

     if ( Pos1 != (int) b1.size( ) && Pos2 != (int) b2.size( ) )
     {    cout << "Warning: strange data encountered in MergeTwoBasevectors.\n";
          PRINT(i);
          PRINT2( pos1, Pos1 );
          PRINT2( pos2, Pos2 );
          PRINT2( b1.size( ), b2.size( ) );
          PrintVisualAlignment( True, cout, b1, b2, a, q1, q2 );    }

     for ( int j = 0; j < (int) b1.size( ) - Pos1; j++ )
     {    if ( i >= static_cast<int>(q.size( )) ) break;         // Don't know if this ever happens.
          if ( j + Pos1 >= static_cast<int>(q1.size( )) ) break; // ditto
          q[i] = q1[j + Pos1];
          c.Set( i++, b1[j + Pos1]);    }
     for ( int j = 0; j < (int) b2.size( ) - Pos2; j++ )
     {    if ( i >= static_cast<int>(q.size( )) ) break;         // Don't know if this ever happens.
          if ( j + Pos2 >= static_cast<int>(q2.size( )) ) break; // ditto;
          q[i] = q2[j + Pos2];
          c.Set( i++, b2[j + Pos2]);    }
     c.resize(i);
     q.resize(i);    }

// AlignTwoBasevectors( b1, b2, a, min_overlap, max_overlap, max_error_rate, log ):
//
// attempt to align basevector b1 to basevector b2, subject to the following
// criteria:
//
// * between min_overlap and max_overlap bases at the right end of b1 must align
//   with the left end of b2;
//
// * the error rate on this overlap must be <= max_error_rate;
//
// * subject to these two conditions, there must be a unique solution.
//
// Return the number of overlapping bases if successful, else -1.
// Also set a to the alignment of the two basevectors.

// Note.  The alignment method that is used is very inefficient.

// Note.  Additional arguments have been added.

class augmented_align {

     public:

     align a;
     int RC;
     int length1, length2;
     int id1, id2;
     int pos1, pos2;
     Float score;
     int Pos1, Pos2;

     augmented_align( ) { }

     void SetFrom( const align& ax, int RCx, int length1x, int length2x, int id1x,
          int id2x, int pos1x, int pos2x, Float scorex, int Pos1x, int Pos2x )
     {    a = ax;
          RC = RCx;
          length1 = length1x;
          length2 = length2x;
          id1 = id1x;
          id2 = id2x;
          pos1 = pos1x;
          pos2 = pos2x;
          score = scorex;
          Pos1 = Pos1x;
          Pos2 = Pos2x;    }

};

int AlignTwoBasevectors( const basevector& b1, const basevector& b2, align& a,
     int min_overlap, int max_overlap, float max_error_rate, ostream* log, int& Rc,
     int mode, int K, int stretch, int nstretch, const qualvector& q1,
     const qualvector& q2, float max_score, float max_errors,
     ostream& badlog, Bool avoid_promiscuous_kmers, int max_cliq8, int max_aligns8,
     int max_err8, int local_max_errors, Bool alt_method, int bandwidth,
     int min_mutmer, float ambiguity_threshold, float max_polyscore_percent,
     Bool sw_gap_method, double min_progression_ratio,
     makealigns_method *method_ptr, Bool assert_if_problem,
     vec< vec< vec<int> > >* allowed_offsets, int max_offset_discrep )

{
     START_TIMER( atb1, 2000 );
     ForceAssert( max_aligns8 <= 50000 ); // else have to change als, below
     ForceAssert( K == 8 || K == 12 || K == 16 || K == 24 );
     ForceAssert( b1.size( ) != 0 );
     ForceAssert( b2.size( ) != 0 );

     vec<bvec const*> EE;
     EE.reserve(2);
     EE.push_back(&b1);
     EE.push_back(&b2);

     if ( stretch == 0 )
     {    if ( K == 8 ) stretch = 8;
          else if ( K == 12 ) stretch = 6;
          else if ( K == 16 ) stretch = 4;
          else if ( K == 24 ) stretch = 2;    }

     int makealigns_max_errs = Max( 1000, (int)max_errors );

     static vec<align_plus> answer;
     static vec<int> answer_counts;
     STOP_TIMER(atb1);
     START_TIMER(atb2, 2000);

     makealigns_method *method;

     if ( method_ptr != 0 )
       method = method_ptr;

     else if ( alt_method )
     {
       static makealigns_alt_method *alt_method = new makealigns_alt_method;

       alt_method->SetBandwidth( bandwidth );
       alt_method->SetMaxErrs( (K==8 ? max_err8 : makealigns_max_errs) );

       method = alt_method;
     }
     else if ( sw_gap_method )
     {
       makealigns_sw_gap_method *sw_gap_method = new makealigns_sw_gap_method;

       sw_gap_method->SetMaxErrs(500);
       sw_gap_method->SetEndStretch(stretch);
       sw_gap_method->SetMaxMutmerOffsetDiff(1000);
       sw_gap_method->SetMaxGap(1000);
       sw_gap_method->SetAffinePenalties(True);
       sw_gap_method->SetMinMaxMutmerLength(100);
       sw_gap_method->SetMinProgressionRatio(min_progression_ratio);

       method = sw_gap_method;    }
     else
     {
       static makealigns_orig_method *orig_method = new makealigns_orig_method;

       orig_method->SetMaxBadness( 100 );
       orig_method->SetMaxErrs( (K==8 ? max_err8 : makealigns_max_errs) );
       orig_method->SetLocalMaxErrs( local_max_errors );
       orig_method->SetLocalMaxErrsDone( 0 );
       orig_method->SetStretch( stretch );
       orig_method->SetEndStretch( (K==8 ? Max( 2, stretch ) : stretch) );
       orig_method->SetNStretch( nstretch );
       orig_method->SetCl( 20 );
       orig_method->SetMaxAlignCtorCalls( 1000000 );

       method = orig_method;
     }

     #define FOR_K(_K) \
         MakeAlignsStatic<2, _K, 50>( answer, answer_counts, 100, 100, EE,        \
             to_compare(FIRST_VS_SECOND, 1), method, _K==8 ? max_cliq8 : 1000,   \
             _K==8 ? max_aligns8 : 50000, min_mutmer, avoid_promiscuous_kmers,   \
             allowed_offsets, max_offset_discrep )

     DISPATCH_ON_K(K, FOR_K);

     STOP_TIMER(atb2);

     int aligns_length;
     int id1, id2;
     static vector<augmented_align> nobbits(100);
     int nobbit_count = 0;

     Bool rc = False;
     static align als[50000];

     int answer_ptr = 0, answer_counts_ptr = 0;
     while(1)
     {    if ( answer_counts_ptr == (int) answer_counts.size( ) ) break;
          aligns_length = answer_counts[ answer_counts_ptr++ ];
          for ( int ll = 0; ll < aligns_length; ll++ )
          {    ForceAssert( ll < 50000 );
               align_plus& ap = answer[ answer_ptr++ ];
               als[ll] = ap.a;
               rc = ap.RC;
               id1 = ap.id1;
               id2 = ap.id2;
               int over;
               if ( mode != 3333 )
                    over = EstimatedOverlap( als[ll], *EE[id1], *EE[id2] );
               else over = als[ll].Pos1( ) - als[ll].pos1( );
               if ( over <= 0 || over < min_overlap ) continue;
               int length1 = EE[id1]->size( ), length2 = EE[id2]->size( );
               int pos1 = als[ll].pos1( ), pos2 = als[ll].pos2( );
               int Pos1 = als[ll].Pos1( ), Pos2 = als[ll].Pos2( );

               if ( nobbit_count >= (int) nobbits.size( ) )
                    nobbits.resize( nobbits.size( ) + nobbits.size( )/2 );

               if ( id1 < id2 )
               {    if ( !rc )
                    {    nobbits[nobbit_count++].SetFrom( als[ll], rc,
                              length1, length2, id1, id2, pos1, pos2,
                              0, Pos1, Pos2 );    }
                    else
                    {    static align r;
                         r = als[ll];
                         r.ReverseThis( length1, length2 );
                         nobbits[nobbit_count++].SetFrom( r, rc, length1,
                              length2, id1, id2, r.pos1( ), r.pos2( ), 0,
                              r.Pos1( ), r.Pos2( ) );    }    }
               else
               {    static align a;
                    a = als[ll];
                    a.Flip( );
                    nobbits[nobbit_count++].SetFrom( a, rc, length2,
                         length1, id2, id1, pos2, pos1, 0, Pos2, Pos1 );  }  }  }

     // Some variables for each nobbit: accept, errors, score, over.

     static vec<Bool> accept;
     accept.resize_and_set( nobbit_count, False );
     static vec<int> over;
     over.resize( nobbit_count );
     static vec<Float> score, errors;
     score.resize( nobbit_count );
     errors.resize( nobbit_count );

     // Determine which nobbits are accepted, and for those, define the other
     // variables.

     for ( int i = 0; i < nobbit_count; i++ )
     {    align& al = nobbits[i].a;
          int RC = nobbits[i].RC;
          int pos1 = nobbits[i].pos1, Pos1 = nobbits[i].Pos1;
          int pos2 = nobbits[i].pos2;

          if ( mode == 1 && ( RC || pos2 > 1 ) ) continue;
          if ( mode == 2 )
          {    if ( pos1 > 1 || Pos1 < (int) b1.size( ) - 1 ) continue;    }
          if ( mode == 3 && ( !RC || pos1 > 1 ) ) continue;
          if ( mode == 4 && ( RC || (pos1 > 1 && pos2 > 1) ) ) continue;
          if ( mode == 5 && RC ) continue;
          if ( mode == 6 && ( RC || ( pos2 > 1 && Pos1 < (int) b1.size( ) - 1 ) ) )
               continue;
          if ( mode == 7 && RC != Rc ) continue;
          if ( mode == 8 && ( RC ) ) continue;
          if ( mode == 9 && ( !RC ) ) continue;
          if ( mode == 3333 )
          {    if (RC) continue;
               if ( pos2 > 0 && pos1 < (int) b1.size( ) ) continue;    }

          if ( mode != 3333 ) over[i] = EstimatedOverlap( al, b1, b2 );
          else over[i] = al.Pos1( ) - al.pos1( );
          if ( over[i] > max_overlap )
          {    if ( log != 0 ) *log << "EXPECTED OVERLAP EXCEEDED\n";
               continue;    }

          static align al_flip;
          al_flip = al;
          al_flip.Flip( );

          errors[i] = ActualErrors( RC, b2, b1, al_flip, 1, 1 );

          if ( q1.size( ) > 0 )
               score[i] = ScoreAlignment(RC, al_flip, b2, q2, b1, q1);
          else score[i] = Float(errors[i]) / Float(over[i]);

          // Test for excessive errors.

          if ( Float(errors[i]) / Float(over[i]) > Float(max_error_rate) ) continue;
          if ( q1.size( ) > 0 )
          {    if ( score[i] > Float(max_score) && errors[i] > Float(max_errors) )
                    continue;    }
          if ( max_polyscore_percent < 100.0 )
          {    int pscore = ScoreAlignmentPoly( RC, al_flip, b2, q2, b1, q1 );
               if ( Float(100.0) * Float(pscore) / Float( al.Pos1( ) - al.pos1( ) )
                    > Float(max_polyscore_percent) ) continue;    }

          // Accept the overlap.

          accept[i] = True;    }

     // Find the lowest scoring overlap.  If none, return.

     int best_i = -1;
     Float best_score = 0.0;
     for ( int i = 0; i < nobbit_count; i++ )
     {    if ( accept[i] )
          {    if ( best_i < 0 || score[i] < best_score )
               {    best_i = i;
                    best_score = score[i];    }    }    }
     if ( best_i < 0 ) return -1;

     // Is there an overlap whose score is no worse than ambiguity_threshold times
     // the best_score, and which appears to be a fundamentally different overlap?
     // If so, the overlap is ambiguous -- return.

     for ( int i = 0; i < nobbit_count; i++ )
     {    if ( accept[i] && i != best_i
               && score[i] <= Float(ambiguity_threshold) * best_score )
          {
	       augmented_align& other = nobbits[i];
	       augmented_align& best = nobbits[best_i];

               // We declare the two alignments to be "the same" if there exists
               // (j,k) such that position j on the first sequence aligns to
               // position k on the second sequence, according to both alignments.

               Bool same = False;
               for ( int p1 = other.pos1; p1 < other.Pos1; ++p1 )
               {    if ( p1 >= best.pos1 && p1 < best.Pos1 )
                    {    int other_p2 = CorrelatePositions( other.a, p1 );
                         int best_p2 = CorrelatePositions( best.a, p1 );
                         if ( other_p2 >= 0 && other_p2 == best_p2 )
                         {    same = True;
                              break;    }    }    }

               if ( !same )
               {    if ( log != 0 ) *log << "Multiple alignments found.\n";
                    return -1;    }    }    }

     // Check for proper overlap.

     a = nobbits[best_i].a;
     if ( assert_if_problem && mode != 3333 )
     {    if ( ( a.pos1( ) > 0 && a.pos2( ) > 0 )
               || ( a.Pos1( ) < (int) b1.size( )
                    && a.Pos2( ) < (int) b2.size( ) ) )
          {    cout << "AlignTwoBasevectors: overlap is improper.\n";
               PRINT2( a.pos1( ), a.pos2( ) );
               PRINT2( a.Pos1( ), b1.size( ) );
               PRINT2( a.Pos2( ), b2.size( ) );
               b1.Print( cout, "b1" );
               b2.Print( cout, "b2" );
               ForceAssert( 0 == 1 );    }    }

     // Return the overlap.

     Rc = nobbits[best_i].RC;
     return over[best_i];    }
