///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "Badness.h"
#include "Basevector.h"
#include "math/Functions.h"
#include "PackAlign.h"
#include "PrintAlignment.h"
#include "ShortVector.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "system/System.h"

// ================================================================================
//
// Given a list of nonoverlapping mutmers (in order), create a align
// from them by filling in the gaps with a restricted Smith-Waterman.  You
// have to also provide the reads themselves.
//
// Set errors_found to the number of errors in the alignment.
//
// If there are >= max_errors, abort the construction, returning an 
// align with Nblocks( ) == 0.
// uninitialized alignment.
//
// Does the code check for overflow in the vectors gaps and lengths??
//
// ================================================================================


// helper functions to append and prepend aligns to gap and length vectors

bool AppendAlignment( avector<int> &gaps, int &gaps_end,
		      avector<int> &lengths, int &lengths_end,
		      int &p1, int &p2,
		      alignment &a )
{
  align gap_align( a );
  
  //  PRINT( gap_align.Nblocks() );

  if ( gap_align.Nblocks() == 0 ) return false;
	    
  int gap_pos1 = gap_align.pos1();
  int gap_pos2 = gap_align.pos2();

  // PRINT( gap_pos1 );
  // PRINT( gap_pos2 );

  ForceAssert( ! ( gap_pos1 > 0 && gap_pos2 > 0 ) );

  if ( gap_pos1 > 0 )
  {
    gaps(gaps_end++) = -gap_pos1;
    p1 += gap_pos1;
  }
  else if ( gap_pos2 > 0 )
  {
    gaps(gaps_end++) = gap_pos2;
    p2 += gap_pos2;
  }
  else
  {
    gaps(gaps_end++) = 0;
  }
	
  ForceAssertEq( gap_align.Gaps(0), 0 );

  int new_length = gap_align.Lengths( 0 );
  lengths(lengths_end++) = new_length;
  p1 += new_length;
  p2 += new_length;
  
  ForceAssertGe( new_length, 0 );

  // PRINT( new_length );

  for ( int block_idx = 1; block_idx < gap_align.Nblocks(); ++block_idx )
  {
    int new_gap = gap_align.Gaps( block_idx );
    gaps(gaps_end++) = new_gap;
    if ( new_gap < 0 )
      p1 += -new_gap;
    else
      p2 += new_gap;
    
    // PRINT( new_gap );
	      
    int new_length = gap_align.Lengths( block_idx );
    lengths(lengths_end++) = new_length;
    p1 += new_length;
    p2 += new_length;
    
    ForceAssertGe( new_length, 0 );

    // PRINT( new_length );
	      
    // PRINT2( p1, p2 );
  }

  return true;
}


bool PrependAlignment( avector<int> &gaps, int &gaps_start,
		       avector<int> &lengths, int &lengths_start,
		       int &pos1, int &pos2,
		       int gap_len1, int gap_len2,
		       alignment &a )
{
  align gap_align( a );

  ForceAssertEq( gap_align.Gaps(0), 0 );
       
  if ( gap_align.Nblocks() == 0 ) return false;
  
  int post_gap1 = gap_len1 - gap_align.Pos1();
  int post_gap2 = gap_len2 - gap_align.Pos2();
  
  ForceAssertGe( post_gap1, 0 );
  ForceAssertGe( post_gap2, 0 );

  // PRINT( post_gap1 );
  // PRINT( post_gap2 );
  
  ForceAssert( post_gap1 == 0 || post_gap2 == 0 );
  
  ForceAssertEq( gaps(gaps_start), 0 );
  
  if ( post_gap1 > 0 )
  {
    gaps(gaps_start) = -post_gap1;
    pos1 -= post_gap1;
  }
  else if ( post_gap2 > 0 )
  {
    gaps(gaps_start) = post_gap2;
    pos2 -= post_gap2;
  }
  else
  {
    gaps(gaps_start) = 0;
  }
  
  for ( int block_ridx = 0; block_ridx < gap_align.Nblocks(); ++block_ridx )
  {
    int block_idx = gap_align.Nblocks() - block_ridx - 1;
    
    int new_length = gap_align.Lengths( block_idx );
    lengths(--lengths_start) = new_length;
    pos1 -= new_length;
    pos2 -= new_length;
    
    ForceAssertGe( new_length, 0 );

    // PRINT( new_length );
    
    int new_gap = gap_align.Gaps( block_idx );
    gaps(--gaps_start) = new_gap;
    if ( new_gap < 0 )
      pos1 -= -new_gap;
    else
      pos2 -= new_gap;
    
    // PRINT2( block_idx, new_gap );
      
    // PRINT2( pos1, pos2 );
  }
  
  ForceAssertEq( gaps(gaps_start), 0 );
  
  return true;
}



void align::CreateFromMutmersAndSW( int k, shortvector<mutmer>& m, 
				    const basevector& rd1, const basevector& rd2, 
				    int max_errors, int end_stretch, int& errors_found,
				    bool affine_penalties )
{    
//   PRINT( m.length );
     SetNblocks(0);
     // cout << "processing " << m.length << " mutmers\n"; // XXX
     int pos1, pos2, len, errors;
     m(0).Unpack(pos1, pos2, len, errors);
     ForceAssertLe( pos1 + len, (int) rd1.size( ) ); // XXX
     ForceAssertLe( pos2 + len, (int) rd2.size( ) ); // XXX
     if ( errors >= max_errors ) 
     {
       errors_found = errors;
       return;
     }

//      PRINT4( pos1, pos2, len, errors );
     
//      PRINT2( rd1.size(), rd2.size() );

     int halfgaps = 10000;
     avector<int> gaps, lengths;
     if ( max_errors > halfgaps )
       halfgaps = max_errors;
     gaps.resize( 2 * halfgaps );
     lengths.resize( 2 * halfgaps );

     int gaps_end = halfgaps, 
          lengths_end = halfgaps;   // gaps_end increases and
     int gaps_start = halfgaps, 
          lengths_start = halfgaps;    // gaps_start decrease in this loop
     gaps(gaps_end++) = 0;
     lengths(lengths_end++) = len;
     int p1 = pos1 + len, p2 = pos2 + len;
     int P1, P2, L, E;
     for ( int i = 1; i < m.length; i++ )
     {    
          m(i).Unpack(P1, P2, L, E);
          ForceAssertLe( P1 + L, (int) rd1.size( ) ); // XXX
          ForceAssertLe( P2 + L, (int) rd2.size( ) ); // XXX
          errors += E;
	  
//  	  PRINT4( P1, P2, L, E );

          if ( errors >= max_errors ) 
	  {
	    errors_found = errors;
	    return;
	  }
          int gap1 = P1 - p1, gap2 = P2 - p2;

//  	  PRINT2( p1, p2 );
//  	  PRINT2( P1, P2 );
//  	  PRINT2( gap1, gap2 );

          // Trim if need be.

          if ( gap1 < 0 )
          {    L += gap1;
               if ( L <= 0 ) continue;

               unsigned int trimmed_P1 = P1 - gap1;
               unsigned int trimmed_P2 = P2 - gap1;
               if ( trimmed_P1 >= rd1.size() || 
                    trimmed_P2 >= rd2.size() ) continue;

               for ( int j = 0; j < -gap1; j++ )
                    if ( rd1[P1 + j] != rd2[P2 + j] ) --errors;

               P1 = trimmed_P1;
               P2 = trimmed_P2;

               gap1 = P1 - p1;
               gap2 = P2 - p2;    }
	  ForceAssertGe( gap1, 0 );

          if ( gap2 < 0 )
          {    L += gap2;
               if ( L <= 0 ) continue;

               unsigned int trimmed_P1 = P1 - gap2;
               unsigned int trimmed_P2 = P2 - gap2;
               if ( trimmed_P1 >= rd1.size() || 
                    trimmed_P2 >= rd2.size() ) continue;

               for ( int j = 0; j < -gap2; j++ )
                    if ( rd1[P1 + j] != rd2[P2 + j] ) --errors;

               P1 = trimmed_P1;
               P2 = trimmed_P2;

               gap1 = P1 - p1;
               gap2 = P2 - p2;    }
	  ForceAssertGe( gap2, 0 );

// 	  PRINT2( p1, p2 );
// 	  PRINT2( gap1, gap2 );

          // Here's the situation now.  Starting at p1 on read 1, and
          // p2 on read 2, we've got gap1 bases on read 1 and gap2
          // bases on read 2.  If one or both are zero, insert an
          // appropriate gap and then the mutmer.  Otherwise, generate
          // alignment data between the two gaps and insert that and
          // then the mutmer.

	  if ( gap1 == 0 && gap2 == 0 )
	  {
	    lengths(lengths_end-1) += L;
	    p1 += L;
	    p2 += L;
	  }
	  else if ( gap1 == 0 )
	  {
	    gaps(gaps_end++) = gap2;
	    p2 += gap2;

	    lengths(lengths_end++) = L;
	    p1 += L;
	    p2 += L;
	  }
	  else if ( gap2 == 0 )
	  {
	    gaps(gaps_end++) = -gap1;
	    p1 += gap1;

	    lengths(lengths_end++) = L;
	    p1 += L;
	    p2 += L;
	  }
	  else
	  {
	    // Expand the region we're interested in by as much as k / 2 - 1
	    // bases on each side of the gap.
	    
	    int gap_edge = Min( k, L, lengths(lengths_end-1) ) / 2 - 1;
	    
	    // PRINT( gap_edge );
	    
	    if ( gap_edge > 0 )
	    {
	      lengths(lengths_end-1) -= gap_edge;
	      p1 -= gap_edge;
	      p2 -= gap_edge;
	      
	      L -= gap_edge;
	      P1 += gap_edge;
	      P2 += gap_edge;
	      
	      gap1 += 2 * gap_edge;
	      gap2 += 2 * gap_edge;
	    }
	    
// 	    PRINT2( gap1, gap2 );
// 	    PRINT2( p1, p2 );
// 	    PRINT2( P1, P2 );
// 	    PRINT( L );

	    // We perform a Smith-Waterman alignment over the gap region
	    // and paste that alignment between the two mutmers.
	    
	    alignment gap_alignment;
	    
	    basevector rd1_gap, rd2_gap;
	    rd1_gap.SetToSubOf( rd1, p1, gap1 );
	    rd2_gap.SetToSubOf( rd2, p2, gap2 );
	    
	    int offset = 0;

// 	    cout << "Performing SW"
// 		 << ": rd1=" << rd1.size()
// 		 << ", rd2=" << rd2.size()
// 		 << ", gap1=" << gap1 
// 		 << ", gap2=" << gap2 
// 		 << ", mutmer_len1=" << lengths(lengths_end-1) 
// 		 << ", mutmer_len2=" << L 
// 		 << " " << flush;
	    
// 	    cout << "Performing SW:" << endl;
// 	    rd1_gap.Print( cout, "rd1_gap" );
// 	    rd2_gap.Print( cout, "rd2_gap" );

	    int gap_errors;

	    if( affine_penalties ) 
	    {
	      if ( gap2 >= gap1 )
	        gap_errors = SmithWatAffine( rd1_gap, rd2_gap, gap_alignment );
	      else
	      {
	        gap_errors = SmithWatAffine( rd2_gap, rd1_gap, gap_alignment );
	        gap_alignment.Flip();
	      }
	    }
	    else 
	    {
	      if ( gap2 >= gap1 )
	        gap_errors = SmithWatFree( rd1_gap, rd2_gap, offset, gap_alignment, true, true );
	      else
	      {
	        gap_errors = SmithWatFree( rd2_gap, rd1_gap, offset, gap_alignment, true, true );
	        gap_alignment.Flip();
	      }
	    }
	    
// 	    PRINT( gap_errors );

// 	    PrintVisualAlignment( False, cout, rd1_gap, rd2_gap, gap_alignment );

	    if ( ! AppendAlignment( gaps, gaps_end,
				    lengths, lengths_end,
				    p1, p2, gap_alignment ) )
	    {
	      errors_found = errors + gap_errors;
	      return;
	    }

	    ForceAssert( ! ( p1 < P1 && p2 < P2 ) );
	    
	    if ( (int)p1 < P1 )
	    {
	      int new_gap = P1 - p1;
	      gaps(gaps_end++) = -new_gap;
	      p1 += new_gap;
	      // PRINT( -new_gap );
	    }
	    else if ( (int)p2 < P2 )
	    {
	      int new_gap = P2 - p2;
	      gaps(gaps_end++) = new_gap;
	      p2 += new_gap;
	      // PRINT( new_gap );
	    }
	    else
	    {
	      int new_gap = 0;
	      gaps(gaps_end++) = new_gap;
	      // PRINT( new_gap );
	    }
	    
	    vector<int> mgg = gap_alignment.MutationsGap1Gap2( rd1_gap, rd2_gap );
	    
	    errors += Sum(mgg);
	    
	    // Add the mutmer.
	    
	    // PRINT( L );
	    
	    lengths(lengths_end++) = L;
	    p1 += L;
	    p2 += L;
	  }	    
	   
	  // PRINT2( p1, p2 );
	  
	  ForceAssertLe( p1, (int) rd1.size( ) ); // XXX
	  ForceAssertLe( p2, (int) rd2.size( ) ); // XXX
     }

     ForceAssertLe( p1, (int) rd1.size( ) ); // XXX
     ForceAssertLe( p2, (int) rd2.size( ) ); // XXX

     // Did we need to force a proper alignment?

     bool forced_proper = false;

     // Now extend the alignment to the right.
     int rspace, rspace1, rspace2;

     rspace1 = rd1.size() - p1;
     rspace2 = rd2.size() - p2;
     rspace = Min( rspace1, rspace2 );
     rspace1 = Min( rspace * 2, rspace1 );
     rspace2 = Min( rspace * 2, rspace2 );

     if ( rspace < end_stretch * k && rspace > 0 )
     {
       alignment gap_alignment;
       
       // Expand the region we're interested in by as much as k / 2 - 1
       // bases.
       
       int gap_edge = Min( k, lengths(lengths_end-1) ) / 2 - 1;
       
       // PRINT( gap_edge );

       if ( gap_edge > 0 )
       {
	 lengths(lengths_end-1) -= gap_edge;
	 p1 -= gap_edge;
	 p2 -= gap_edge;
	 
	 rspace1 += gap_edge;
	 rspace2 += gap_edge;
       }

       basevector rd1_gap, rd2_gap;
       rd1_gap.SetToSubOf( rd1, p1, rspace1 );
       rd2_gap.SetToSubOf( rd2, p2, rspace2 );
	    
       int offset = 0;
       
//        cout << "Performing SW to extend to right"
// 	    << ": rd1=" << rd1.size()
// 	    << ", rd2=" << rd2.size()
// 	    << ", gap1=" << rd1_gap.size()
// 	    << ", gap2=" << rd2_gap.size()
// 	    << " " << flush;

       // SmithWatAffine does not support unrestricted edges, so we use SmithWatFree 
       // here for now.
       if ( rspace2 >= rspace1 )
	 SmithWatFree( rd1_gap, rd2_gap, offset, gap_alignment, true, false );
       else
       {
	 SmithWatFree( rd2_gap, rd1_gap, offset, gap_alignment, true, false );
	 gap_alignment.Flip();
       }

//        int gap_errors = gap_alignment.Errors();
	    
// 	    cout << gap_errors << " errors." << endl;

       AppendAlignment( gaps, gaps_end,
			lengths, lengths_end,
			p1, p2, gap_alignment );

       vector<int> mgg = gap_alignment.MutationsGap1Gap2( rd1_gap, rd2_gap );
       
       errors += Sum(mgg);

       // Force it to be proper.

       int basesleft = Min( rd1.size() - p1, rd2.size() - p2 );
       
       if ( basesleft > 0 )
       {
         cout << "Forcing proper alignment of " << basesleft << "bp hanging end on right." << endl;
         PRINT3( rspace1, rspace2, gap_edge );
         lengths( lengths_end-1 ) += basesleft;
         p1 += basesleft;
         p2 += basesleft;
         forced_proper = true;
       }
     }

     //  	  PRINT2( p1, p2 );
     
     ForceAssertLe( p1, (int) rd1.size( ) ); // XXX
     ForceAssertLe( p2, (int) rd2.size( ) ); // XXX

     // Now extend the alignment to the left.

     int lspace, lspace1, lspace2;

     lspace1 = pos1;
     lspace2 = pos2;
     lspace = Min( lspace1, lspace2 );
     lspace1 = Min( lspace * 2, lspace1 );
     lspace2 = Min( lspace * 2, lspace2 );

     if ( lspace < end_stretch * k && lspace > 0 )
     {
       alignment gap_alignment;

       // Expand the region we're interested in by as much as k / 2 - 1
       // bases.
       
       int gap_edge = Min( k, lengths(lengths_start) ) / 2 - 1;
       
       // PRINT( gap_edge );

       if ( gap_edge > 0 )
       {
	 lengths(lengths_start) -= gap_edge;
	 pos1 += gap_edge;
	 pos2 += gap_edge;
	 
	 lspace1 += gap_edge;
	 lspace2 += gap_edge;
       }

       basevector rd1_gap, rd2_gap;
       rd1_gap.SetToSubOf( rd1, pos1 - lspace1, lspace1 );
       rd2_gap.SetToSubOf( rd2, pos2 - lspace2, lspace2 );
	    
       int offset = 0;
       
//        cout << "Performing SW to extend to left"
// 	    << ": rd1=" << rd1.size()
// 	    << ", rd2=" << rd2.size()
//  	    << ", gap1=" << rd1_gap.ToString()
//  	    << ", gap2=" << rd2_gap.ToString() << endl;
// 	    << " " << flush;

       // SmithWatAffine does not support unrestricted edges, so we use SmithWatFree 
       // here for now.	    
       if ( lspace2 >= lspace1 )
	 SmithWatFree( rd1_gap, rd2_gap, offset, gap_alignment, false, true );
       else
       {
	 SmithWatFree( rd2_gap, rd1_gap, offset, gap_alignment, false, true );
	 gap_alignment.Flip();
       }

//        int gap_errors = gap_alignment.Errors();
//        cout << gap_errors << " errors." << endl;

       PrependAlignment( gaps, gaps_start,
			 lengths, lengths_start,
			 pos1, pos2, 
			 lspace1, lspace2, 
			 gap_alignment );

       vector<int> mgg = gap_alignment.MutationsGap1Gap2( rd1_gap, rd2_gap );
       
       errors += Sum(mgg);

       // Force it to be proper.

       int basesleft = Min( pos1, pos2 );
       
       if ( basesleft > 0 )
       {
         cout << "Forcing proper alignment of " << basesleft << "bp hanging end on left." << endl;
         PRINT3( lspace1, lspace2, gap_edge );
         lengths( lengths_start ) += basesleft;
         pos1 -= basesleft;
         pos2 -= basesleft;
         forced_proper = true;
       }
     }

     ForceAssertGe( pos1, 0 ); // XXX
     ForceAssertGe( pos2, 0 ); // XXX
     
     ForceAssertEq( gaps_end - gaps_start, lengths_end - lengths_start );

     ForceAssertLe( p1, (int) rd1.size( ) ); // XXX
     ForceAssertLe( p2, (int) rd2.size( ) ); // XXX

     ForceAssertGt( gaps_start, 1 ) ; 
     ForceAssertLt( gaps_end, 2 * halfgaps - 1 ) ; 

     avector<int> gaps2(2 * halfgaps), lengths2(2 * halfgaps);
     int l = 0;
     for ( int j = gaps_start; j < gaps_end; j++ )
     {
//        cout << "gaps(" << j << ") = " << gaps(j) << endl;
//        cout << "lengths(" << j << ") = " << lengths(j) << endl;

       if ( j == gaps_start || gaps(j) != 0 )
       {
	 ForceAssertGe( lengths(j), 0 );
	 gaps2(l) = gaps(j);
	 lengths2(l) = lengths(j);
	 ++l;
       }
       else 
	 lengths2( l - 1 ) += lengths(j);    
     }

//       PRINT( pos1 );
//       PRINT( pos2 );

     Setpos1(pos1);
     Setpos2(pos2);
     SetNblocks(l);
     for ( int i = 0; i < l; i++ )
     {
       SetGap( i, gaps2(i) );
       SetLength( i, lengths2(i) );    
     }
     errors_found = errors;    


     if ( forced_proper )
     {
       PRINT2( pos1, rd1.size() - Pos1() );
       PRINT2( pos2, rd2.size() - Pos2() );
       PrintVisualAlignment( False, cout, rd1, rd2, *this );
     }
}
