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
#include "ShortVector.h"
#include "system/System.h"

// ================================================================================
//
// Given a list of nonoverlapping mutmers (in order), create a align from
// them by filling in the gaps.  You have to also provide the reads themselves.
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

void align::CreateFromMutmers(int k, shortvector<mutmer>& m, const basevector& rd1, 
     const basevector& rd2, int max_errors, float max_badness,
     int local_max_errors, int end_stretch, int local_max_errors_done, 
     int& errors_found )
{    
     SetNblocks(0);
     const int MAX_ERRORS = 10000 ;
     if ( max_errors > MAX_ERRORS )
          FatalErr( "You've called the align constructor (in "
               << "AlignFromMutmers.cc) with a max_errors\nvalue of > 10000.\n"
               << "If you really want to do this, you need to increase  "
		    << "MAX_ERRORS and halfgaps in AlignFromMutmers.cc\n " );
     // cout << "processing " << m.length << " mutmers\n"; // XXX
     int pos1, pos2, len, errors;
     m(0).Unpack(pos1, pos2, len, errors);
     Assert( pos1 + len <= (int) rd1.size( ) ); // XXX
     Assert( pos2 + len <= (int) rd2.size( ) ); // XXX
     if ( errors >= max_errors ) return;

     #ifdef KEN
          const int halfgaps = 500000 ; 
     #else
          const int halfgaps = MAX_ERRORS ;  // I suspect that this works but I 
                                             // don't want to introduce bugs 
                                             // arbitrarily.
     #endif

     avector<int> gaps(2 * halfgaps), lengths(2 * halfgaps);
     int gaps_length = halfgaps, 
          lengths_length = halfgaps;   // gaps_length increases and
     int gaps_start = halfgaps, 
          lengths_start = halfgaps;    // gaps_start decrease in this loop
     gaps(gaps_length++) = 0;
     lengths(lengths_length++) = len;
     unsigned int p1 = pos1 + len, p2 = pos2 + len;
     int P1, P2, L, E;
     for ( int i = 1; i < m.length; i++ )
     {    m(i).Unpack(P1, P2, L, E);
          Assert( P1 + L <= (int) rd1.size( ) ); // XXX
          Assert( P2 + L <= (int) rd2.size( ) ); // XXX
          errors += E;
          if ( errors >= max_errors ) return;
          int gap1 = P1 - p1, gap2 = P2 - p2;

          // Trim if need be.

          if ( gap1 < 0 )
          {    for ( int j = 0; j < -gap1; j++ )
                    if ( rd1[P1 + j] != rd2[P2 + j] ) --errors;
               P1 -= gap1;
               P2 -= gap1;
               L += gap1;
               gap1 = P1 - p1;
               gap2 = P2 - p2;    }
          if ( gap2 < 0 )
          {    for ( int j = 0; j < -gap2; j++ )
                    if ( rd1[P1 + j] != rd2[P2 + j] ) --errors;
               P1 -= gap2;
               P2 -= gap2;
               L += gap2;
               gap1 = P1 - p1;
               gap2 = P2 - p2;    }

          // Here's the situation now.  Starting at p1 on read 1, and p2 on
          // read 2, we've got gap1 bases on read 1 and gap2 bases on read 2.
          // Generate alignment data for these.

          int errors_was = errors;

          while(1)
          {    int mtch, mtch1, mtch2;

               // Compute the number of matching characters after the first 
               // mismatch.  Allow one extra mismatch.

               int mismtch = 0;
               int first_bug = 0;
               // for ( mtch = 1; mtch < Min(gap1, gap2); mtch++ ) 
               for ( mtch = 1; 
                    mtch < int(Min(rd1.size( ) - p1, rd2.size( ) - p2)); mtch++ )
                    if ( rd1[p1+mtch] != rd2[p2+mtch] )
                    {    if ( mismtch == 1 ) break;
                         first_bug = mtch - 1;
                         ++mismtch;    }
               --mtch;

               // Compute the number of matching characters if we assume a gap 
               // on rd1.  Allow one extra mismatch.

               int mismtch1 = 0;
               for ( mtch1 = 0; mtch1 + 1 < gap1 && mtch1 < gap2; mtch1++ )
                    if ( rd1[p1+mtch1+1] != rd2[p2+mtch1] )
                    {    if ( mismtch1 == 1 ) break;
                         ++mismtch1;    }
               if ( mtch1 > 0 && mismtch1 > 0 
                    && rd1[p1+mtch1] != rd2[p2+mtch1-1] )
               {    --mtch1;
                    --mismtch1;    }

               // Compute the number of matching characters if we assume a 
               // double gap on rd1.  Allow one extra mismatch.

               int misdmtch1 = 0, dmtch1;
               for ( dmtch1 = 0; dmtch1 + 2 < gap1 && dmtch1 < gap2; dmtch1++ )
                    if ( rd1[p1+dmtch1+2] != rd2[p2+dmtch1] )
                    {    if ( misdmtch1 == 1 ) break;
                         ++misdmtch1;    }
               if ( dmtch1 > 0 && misdmtch1 > 0 
                    && rd1[p1+dmtch1+1] != rd2[p2+dmtch1-1] )
               {    --dmtch1;
                    --misdmtch1;    }

               // Compute the number of matching characters if we assume a gap
               // on read2.  Allow one extra mismatch.

               int mismtch2 = 0;
               for ( mtch2 = 0; mtch2 < gap1 && mtch2 + 1 < gap2; mtch2++ )
                    if ( rd1[p1+mtch2] != rd2[p2+mtch2+1] )
                    {    if ( mismtch2 == 1 ) break;
                         ++mismtch2;    }
               if ( mtch2 > 0 && mismtch2 > 0 
                    && rd1[p1+mtch2-1] != rd2[p2+mtch2] )
               {    --mtch2;
                    --mismtch2;    }

               // Compute the number of matching characters if we assume a
               // double gap on read2.  Allow one extra mismatch.

               int misdmtch2 = 0, dmtch2;
               for ( dmtch2 = 0; dmtch2 < gap1 && dmtch2 + 2 < gap2; dmtch2++ )
                    if ( rd1[p1+dmtch2] != rd2[p2+dmtch2+2] )
                    {    if ( misdmtch2 == 1 ) break;
                         ++misdmtch2;    }
               if ( dmtch2 > 0 && misdmtch2 > 0 
                    && rd1[p1+dmtch2-1] != rd2[p2+dmtch2+1] )
               {    --dmtch2;
                    --misdmtch2;    }
               int mtchx = Min( mtch, Min(gap1, gap2) - 1 );
               if ( (mtch - mismtch > 0 || Min(gap1, gap2) == 1) 
                    && mtchx >= 0 && mtch - mismtch >= mtch1 - mismtch1 - 1 
                    && mtch - mismtch >= mtch2 - mismtch2 - 1
                    && mtch - mismtch >= dmtch1 - misdmtch1 - 3
                    && mtch - mismtch >= dmtch2 - misdmtch2 - 3 )
               {    if ( mtchx <= first_bug ) mismtch = 0;
                    mtch = mtchx;
                    lengths( lengths_length - 1 ) += mtch + 1;
                    p1 += mtch + 1;
                    p2 += mtch + 1;
                    gap1 -= mtch + 1;
                    gap2 -= mtch + 1;
                    errors += 1 + mismtch;    }
               else if ( mtch1 - mismtch1 > 0 
                    && mtch1 - mismtch1 >= mtch2 - mismtch2
                    && mtch1 - mismtch1 >= dmtch1 - misdmtch1 - 2
                    && mtch1 - mismtch1 >= dmtch2 - misdmtch2 - 2 )
               {    lengths(lengths_length++) = mtch1;
                    gaps(gaps_length++) = -1;
                    p1 += mtch1 + 1;
                    p2 += mtch1;
                    gap1 -= mtch1 + 1;
                    gap2 -= mtch1;
                    errors += 2 + mismtch1;    }
               else if ( mtch2 - mismtch2 > 0
                    && mtch2 - mismtch2 >= dmtch1 - misdmtch1 - 2
                    && mtch2 - mismtch2 >= dmtch2 - misdmtch2 - 2 )
               {    lengths(lengths_length++) = mtch2;
                    gaps(gaps_length++) = 1;
                    p1 += mtch2;
                    p2 += mtch2 + 1;
                    gap1 -= mtch2;
                    gap2 -= mtch2 + 1;
                    errors += 2 + mismtch2;    }
               else if ( dmtch1 - misdmtch1 > 2
                    && dmtch1 - misdmtch1 >= dmtch2 - misdmtch2 )
               {    lengths(lengths_length++) = dmtch1;
                    gaps(gaps_length++) = -2;
                    p1 += dmtch1 + 2;
                    p2 += dmtch1;
                    gap1 -= dmtch1 + 2;
                    gap2 -= dmtch1;
                    errors += 4 + misdmtch1;    }
               else if ( dmtch2 - misdmtch2 > 2 )
               {    lengths(lengths_length++) = dmtch2;
                    gaps(gaps_length++) = 2;
                    p1 += dmtch2;
                    p2 += dmtch2 + 2;
                    gap1 -= dmtch2;
                    gap2 -= dmtch2 + 2;
                    errors += 4 + misdmtch2;    }
               else if ( gap2 == 0 )
               {    p1 += gap1;
                    while ( gap1 > 19 )
                    {    gaps(gaps_length++) = -19;
                         lengths(lengths_length++) = 0;
                         gap1 -= 19;
                         errors += 38;    }
                    gaps(gaps_length++) = -gap1;
                    lengths(lengths_length++) = L;
                    errors += 2 * gap1;    
                    break;    }
               else if ( gap1 == 0 )
               {    p2 += gap2;
                    while ( gap2 > 19 )
                    {    gaps(gaps_length++) = 19;
                         lengths(lengths_length++) = 0;
                         gap2 -= 19;
                         errors += 38;    }
                    gaps(gaps_length++) = gap2;
                    lengths(lengths_length++) = L;
                    errors += 2 * gap2;    
                    break;    }
               else
               {    --gap1;
                    --gap2;
                    ++p1;
                    ++p2;
                    ++lengths(lengths_length-1);
                    ++errors;    }    }

          p1 += L;
          p2 += L;    
          if ( errors - errors_was > local_max_errors ) return;
          Assert( p1 <= rd1.size( ) ); // XXX
          Assert( p2 <= rd2.size( ) ); // XXX
          }

     // Now walk over any single base errors on the right.
     Assert( p1 <= rd1.size( ) ); // XXX
     Assert( p2 <= rd2.size( ) ); // XXX

     unsigned int rspace;
     int lspace;
     int errors_was = errors;

     rwalk:
     if ( errors - errors_was > local_max_errors_done ) goto done_walking;
     if ( errors - errors_was > local_max_errors ) return;
     if ( errors >= max_errors ) goto done_walking;
     rspace = Min( rd1.size( ) - p1, rd2.size( ) - p2 );
     if ( rspace >= 1 )
     {    int mtch, mtch1, mtch2;

          // Compute the number of matching characters after the first mismatch.
          // Allow one extra mismatch.
          int mismtch = 0;

          for ( mtch = 1; mtch < int(rspace); mtch++ ) 
               if ( rd1[p1+mtch] != rd2[p2+mtch] )
               {    if ( mismtch == 1 ) break;
                    ++mismtch;    }
          --mtch;

          // Compute the number of matching characters if we assume a gap on
          // rd1.  Allow one extra mismatch.
          int mismtch1 = 0;
          for ( mtch1 = 0; 
               p1 + mtch1 + 1 < rd1.size( ) && p2 + mtch1 < rd2.size( ); 
               mtch1++ )
               if ( rd1[p1+mtch1+1] != rd2[p2+mtch1] )
               {    if ( mismtch1 == 1 ) break;
                    ++mismtch1;    }
          if ( mtch1 > 0 && mismtch1 > 0 
               && rd1[p1+mtch1] != rd2[p2+mtch1-1] )
          {    --mtch1;
               --mismtch1;    }

          // Compute the number of matching characters if we assume a double 
          // gap on rd1.  Allow one extra mismatch.
          int misdmtch1 = 0, dmtch1;
          for ( dmtch1 = 0; p1 + dmtch1 + 2 < rd1.size( ) 
               && p2 + dmtch1 < rd2.size( ); dmtch1++ )
               if ( rd1[p1+dmtch1+2] != rd2[p2+dmtch1] )
               {    if ( misdmtch1 == 1 ) break;
                    ++misdmtch1;    }
          if ( dmtch1 > 0 && misdmtch1 > 0 
               && rd1[p1+dmtch1+1] != rd2[p2+dmtch1-1] )
          {    --dmtch1;
               --misdmtch1;    }

          // Compute the number of matching characters if we assume a gap on
          // read2.  Allow one extra mismatch.
          int mismtch2 = 0;
          for ( mtch2 = 0; 
               p1 + mtch2 < rd1.size( ) && p2 + mtch2 + 1 < rd2.size( ); 
               mtch2++ )
               if ( rd1[p1+mtch2] != rd2[p2+mtch2+1] )
               {    if ( mismtch2 == 1 ) break;
                    ++mismtch2;    }
          if ( mtch2 > 0 && mismtch2 > 0 
               && rd1[p1+mtch2-1] != rd2[p2+mtch2] )
          {    --mtch2;
               --mismtch2;    }

          // Compute the number of matching characters if we assume a double
          // gap on read2.  Allow one extra mismatch.
          int misdmtch2 = 0, dmtch2;
          for ( dmtch2 = 0; p1 + dmtch2 < rd1.size( ) 
               && p2 + dmtch2 + 2 < rd2.size( ); dmtch2++ )
               if ( rd1[p1+dmtch2] != rd2[p2+dmtch2+2] )
               {    if ( misdmtch2 == 1 ) break;
                    ++misdmtch2;    }
          if ( dmtch2 > 0 && misdmtch2 > 0 
               && rd1[p1+dmtch2-1] != rd2[p2+dmtch2+1] )
          {    --dmtch2;
               --misdmtch2;    }

          if ( (mtch - mismtch > 0 || rspace == 1) 
               && mtch - mismtch >= mtch1 - mismtch1 - 1 
               && mtch - mismtch >= mtch2 - mismtch2 - 1
               && mtch - mismtch >= dmtch1 - misdmtch1 - 3
               && mtch - mismtch >= dmtch2 - misdmtch2 - 3 )
          {    lengths( lengths_length - 1 ) += mtch + 1;
               p1 += mtch + 1;
               p2 += mtch + 1;
               errors += 1 + mismtch;
               goto rwalk;    }
          if ( mtch1 - mismtch1 > 0 
               && (mtch1 - mismtch1 > mtch2 - mismtch2 ||
                    (mtch1 - mismtch1 == mtch2 - mismtch2 && mismtch1 <= mismtch2))
               && mtch1 - mismtch1 >= dmtch1 - misdmtch1 - 2
               && mtch1 - mismtch1 >= dmtch2 - misdmtch2 - 2 )
          {    lengths(lengths_length++) = mtch1;
               gaps(gaps_length++) = -1;
               p1 += mtch1 + 1;
               p2 += mtch1;
               errors += 2 + mismtch1;
               goto rwalk;    }
          if ( mtch2 - mismtch2 > 0 && mtch2 - mismtch2 >= dmtch1 - misdmtch1 - 2
               && mtch2 - mismtch2 >= dmtch2 - misdmtch2 - 2 )
          {    lengths(lengths_length++) = mtch2;
               gaps(gaps_length++) = 1;
               p1 += mtch2;
               p2 += mtch2 + 1;
               errors += 2 + mismtch2;
               goto rwalk;    }
          if ( dmtch1 - misdmtch1 > 2 && dmtch1 - misdmtch1 >= dmtch2 - misdmtch2 )
          {    lengths(lengths_length++) = dmtch1;
               gaps(gaps_length++) = -2;
               p1 += dmtch1 + 2;
               p2 += dmtch1;
               errors += 4 + misdmtch1;
               goto rwalk;    }
          if ( dmtch2 - misdmtch2 > 2 )
          {    lengths(lengths_length++) = dmtch2;
               gaps(gaps_length++) = 2;
               p1 += dmtch2;
               p2 += dmtch2 + 2;
               errors += 4 + misdmtch2;
               goto rwalk;    }
          if ( int(rspace) < end_stretch*k )
          {    ++p1;
               ++p2;
               ++lengths(lengths_length-1);
               ++errors;    
               goto rwalk;    }    }

     // Now walk over any single base errors on the left.
     Assert( p1 <= rd1.size( ) ); // XXX
     Assert( p2 <= rd2.size( ) ); // XXX
     errors_was = errors;

     lwalk:
     // PRINT(errors); // XXX
     if ( errors - errors_was > local_max_errors_done ) goto done_walking;
     if ( errors - errors_was > local_max_errors ) return;
     if ( errors >= max_errors ) goto done_walking;
     lspace = Min( pos1, pos2 );
     if ( lspace >= 1 )
     {    int mtch, mtch1, dmtch1, mtch2, dmtch2;

          // Compute the number of matching characters before the first 
          // mismatch.  Allow one extra mismatch.
          int mismtch = 0;
          int mtchx = 0;
          for ( mtch = 1; mtch + 1 <= lspace; mtch++ ) 
               if ( rd1[pos1-mtch-1] != rd2[pos2-mtch-1] )
               {    if ( mismtch == 1 ) break;
                    mtchx = mtch - 1;
                    ++mismtch;    }
          --mtch;

          // Compute the number of matching characters if we assume a gap on
          // rd1.  Allow one extra mismatch.
          int mismtch1 = 0, mtch1x = 0;
          for ( mtch1 = 1; pos1 - mtch1 >= 0 && pos2-mtch1-1 >= 0; mtch1++ )
               if ( rd1[pos1-mtch1] != rd2[pos2-mtch1-1] )
               {    if ( mismtch1 == 1 ) break;
                    mtch1x = mtch1 - 1;
                    ++mismtch1;    }
          --mtch1;

          // Compute the number of matching characters if we assume a double 
          // gap on rd1.  Allow one extra mismatch.
          int misdmtch1 = 0, dmtch1x = 0;
          for ( dmtch1 = 1; pos1 - dmtch1 >= 0 && pos2-dmtch1-2 >= 0; dmtch1++ )
               if ( rd1[pos1-dmtch1] != rd2[pos2-dmtch1-2] )
               {    if ( misdmtch1 == 1 ) break;
                    dmtch1x = dmtch1 - 1;
                    ++misdmtch1;    }
          --dmtch1;

          // Compute the number of matching characters if we assume a gap on
          // read2.  Allow one extra mismatch.
          int mismtch2 = 0, mtch2x = 0;
          for ( mtch2 = 1; pos1 - mtch2 - 1 >= 0 && pos2 - mtch2 >= 0; mtch2++ )
               if ( rd1[pos1-mtch2-1] != rd2[pos2-mtch2] )
               {    if ( mismtch2 == 1 ) break;
                    mtch2x = mtch2 - 1;
                    ++mismtch2;    }
          --mtch2;

          // Compute the number of matching characters if we assume a double 
          // gap on read2.
          int misdmtch2 = 0, dmtch2x = 0;
          for ( dmtch2 = 1; pos1-dmtch2-2 >= 0 && pos2 - dmtch2 >= 0; dmtch2++ )
               if ( rd1[pos1-dmtch2-2] != rd2[pos2-dmtch2] )
               {    if ( misdmtch2 == 1 ) break;
                    dmtch2x = dmtch2 - 1;
                    ++misdmtch2;    }
          --dmtch2;

          if ( (mtch - mismtch > 0 || lspace == 1) 
               && mtch - mismtch >= mtch1 - mismtch1 - 1 
                    && mtch - mismtch >= mtch2 - mismtch2 - 1 
               && mtch - mismtch >= dmtch1 - misdmtch1 - 3 
               && mtch - mismtch >= dmtch2 - misdmtch2 - 3 )
          {    if ( mismtch == 1 && mtchx > 0 ) 
               {    mtch = mtchx;
                    mismtch = 0;    }
               lengths(lengths_start) += mtch + 1;
               pos1 -= mtch + 1;
               pos2 -= mtch + 1;
               errors += 1 + mismtch;
               goto lwalk;    }
          if ( mtch1 - mismtch1 > 0 && mtch1 - mismtch1 >= mtch2 - mismtch2 
               && mtch1 - mismtch1 >= dmtch1 - misdmtch1 - 2 &&
               mtch1 - mismtch1 >= dmtch2 - misdmtch2 - 2 )
          {    if ( mismtch1 == 1 && mtch1x > 0 ) 
               {    mtch1 = mtch1x;
                    mismtch1 = 0;    }
               lengths(--lengths_start) = mtch1;
               gaps(gaps_start) = 1;
               gaps(--gaps_start) = 0;
               pos1 -= mtch1;
               pos2 -= mtch1 + 1;
               errors += 2 + mismtch1;
               goto lwalk;    }
          if ( mtch2 - mismtch2 > 0 
               && mtch2 - mismtch2 >= dmtch1 - misdmtch1 - 2 
               && mtch2 - mismtch2 >= dmtch2 - misdmtch2 - 2 )
          {    if ( mismtch2 == 1 && mtch2x > 0 ) 
               {    mtch2 = mtch2x;
                    mismtch2 = 0;    }
               lengths(--lengths_start) = mtch2;
               gaps(gaps_start) = -1;
               gaps(--gaps_start) = 0;
               pos1 -= mtch2 + 1;
               pos2 -= mtch2;
               errors += 2 + mismtch2;
               goto lwalk;    }    
          if ( dmtch1 - misdmtch1 > 2 
               && dmtch1 - misdmtch1 >= dmtch2 - misdmtch2 )
          {    if ( misdmtch1 == 1 && dmtch1x > 0 ) 
               {    dmtch1 = dmtch1x;
                    misdmtch1 = 0;    }
               lengths(--lengths_start) = dmtch1;
               gaps(gaps_start) = 2;
               gaps(--gaps_start) = 0;
               pos1 -= dmtch1;
               pos2 -= dmtch1 + 2;
               errors += 4;
               goto lwalk;    }
          if ( dmtch2 - misdmtch2 > 2 )
          {    if ( misdmtch2 == 1 && dmtch2x > 0 ) 
               {    dmtch2 = dmtch2x;
                    misdmtch2 = 0;    }
               lengths(--lengths_start) = dmtch2;
               gaps(gaps_start) = -2;
               gaps(--gaps_start) = 0;
               pos1 -= dmtch2 + 2;
               pos2 -= dmtch2;
               errors += 4;
               goto lwalk;    }    
          if ( lspace < end_stretch*k )
          {    --pos1;
               --pos2;
               if ( gaps(gaps_start) == 0 ) ++lengths(lengths_start);
               else
               {    gaps(--gaps_start) = 0;
                    lengths(--lengths_start) = 1;    }
               ++errors;
               goto lwalk;    }    }

     done_walking:
     Assert( gaps_length - gaps_start == lengths_length - lengths_start );

     Assert( p1 <= rd1.size( ) ); // XXX
     Assert( p2 <= rd2.size( ) ); // XXX

     Assert( gaps_start > 1 ) ; 
     Assert( gaps_length < 2 * halfgaps - 1 ) ; 

     avector<int> gaps2(2 * halfgaps), lengths2(2 * halfgaps);
     int l = 0;
     for ( int j = gaps_start; j < gaps_length; j++ )
     {    if ( j == gaps_start || gaps(j) != 0 )
          {    gaps2(l) = gaps(j);
               lengths2(l) = lengths(j);
               ++l;    }
          else lengths2( l - 1 ) += lengths(j);    }

     Setpos1(pos1);
     Setpos2(pos2);
     SetNblocks(l);
     for ( int i = 0; i < l; i++ )
     {    SetGap( i, gaps2(i) );
          SetLength( i, lengths2(i) );    }
     errors_found = errors;    }

void align::CreateFromMutmersMT(int k, shortvector<mutmer>& m, const basevector& rd1, 
     const basevector& rd2, int max_errors, float max_badness,
     int local_max_errors, int end_stretch, int local_max_errors_done, 
     int& errors_found )
{    
     SetNblocks(0);
     const int MAX_ERRORS = 10000 ;
     if ( max_errors > MAX_ERRORS )
          FatalErr( "You've called the align constructor (in "
               << "AlignFromMutmers.cc) with a max_errors\nvalue of > 10000.\n"
               << "If you really want to do this, you need to increase  "
		    << "MAX_ERRORS and halfgaps in AlignFromMutmers.cc\n " );
     // cout << "processing " << m.length << " mutmers\n"; // XXX
     int pos1, pos2, len, errors;
     m(0).Unpack(pos1, pos2, len, errors);
     Assert( pos1 + len <= (int) rd1.size( ) ); // XXX
     Assert( pos2 + len <= (int) rd2.size( ) ); // XXX
     if ( errors >= max_errors ) return;

     #ifdef KEN
          const int halfgaps = 500000 ; 
     #else
          const int halfgaps = MAX_ERRORS ;  // I suspect that this works but I 
                                             // don't want to introduce bugs 
                                             // arbitrarily.
     #endif

     avector<int> gaps(2 * halfgaps), lengths(2 * halfgaps);
     int gaps_length = halfgaps, 
          lengths_length = halfgaps;   // gaps_length increases and
     int gaps_start = halfgaps, 
          lengths_start = halfgaps;    // gaps_start decrease in this loop
     gaps(gaps_length++) = 0;
     lengths(lengths_length++) = len;
     unsigned int p1 = pos1 + len, p2 = pos2 + len;
     int P1, P2, L, E;
     for ( int i = 1; i < m.length; i++ )
     {    m(i).Unpack(P1, P2, L, E);
          Assert( P1 + L <= (int) rd1.size( ) ); // XXX
          Assert( P2 + L <= (int) rd2.size( ) ); // XXX
          errors += E;
          if ( errors >= max_errors ) return;
          int gap1 = P1 - p1, gap2 = P2 - p2;

          // Trim if need be.

          if ( gap1 < 0 )
          {    for ( int j = 0; j < -gap1; j++ )
                    if ( rd1[P1 + j] != rd2[P2 + j] ) --errors;
               P1 -= gap1;
               P2 -= gap1;
               L += gap1;
               gap1 = P1 - p1;
               gap2 = P2 - p2;    }
          if ( gap2 < 0 )
          {    for ( int j = 0; j < -gap2; j++ )
                    if ( rd1[P1 + j] != rd2[P2 + j] ) --errors;
               P1 -= gap2;
               P2 -= gap2;
               L += gap2;
               gap1 = P1 - p1;
               gap2 = P2 - p2;    }

          // Here's the situation now.  Starting at p1 on read 1, and p2 on
          // read 2, we've got gap1 bases on read 1 and gap2 bases on read 2.
          // Generate alignment data for these.

          int errors_was = errors;

          while(1)
          {    int mtch, mtch1, mtch2;

               // Compute the number of matching characters after the first 
               // mismatch.  Allow one extra mismatch.

               int mismtch = 0;
               int first_bug = 0;
               // for ( mtch = 1; mtch < Min(gap1, gap2); mtch++ ) 
               for ( mtch = 1; 
                    mtch < int(Min(rd1.size( ) - p1, rd2.size( ) - p2)); mtch++ )
                    if ( rd1[p1+mtch] != rd2[p2+mtch] )
                    {    if ( mismtch == 1 ) break;
                         first_bug = mtch - 1;
                         ++mismtch;    }
               --mtch;

               // Compute the number of matching characters if we assume a gap 
               // on rd1.  Allow one extra mismatch.

               int mismtch1 = 0;
               for ( mtch1 = 0; mtch1 + 1 < gap1 && mtch1 < gap2; mtch1++ )
                    if ( rd1[p1+mtch1+1] != rd2[p2+mtch1] )
                    {    if ( mismtch1 == 1 ) break;
                         ++mismtch1;    }
               if ( mtch1 > 0 && mismtch1 > 0 
                    && rd1[p1+mtch1] != rd2[p2+mtch1-1] )
               {    --mtch1;
                    --mismtch1;    }

               // Compute the number of matching characters if we assume a 
               // double gap on rd1.  Allow one extra mismatch.

               int misdmtch1 = 0, dmtch1;
               for ( dmtch1 = 0; dmtch1 + 2 < gap1 && dmtch1 < gap2; dmtch1++ )
                    if ( rd1[p1+dmtch1+2] != rd2[p2+dmtch1] )
                    {    if ( misdmtch1 == 1 ) break;
                         ++misdmtch1;    }
               if ( dmtch1 > 0 && misdmtch1 > 0 
                    && rd1[p1+dmtch1+1] != rd2[p2+dmtch1-1] )
               {    --dmtch1;
                    --misdmtch1;    }

               // Compute the number of matching characters if we assume a gap
               // on read2.  Allow one extra mismatch.

               int mismtch2 = 0;
               for ( mtch2 = 0; mtch2 < gap1 && mtch2 + 1 < gap2; mtch2++ )
                    if ( rd1[p1+mtch2] != rd2[p2+mtch2+1] )
                    {    if ( mismtch2 == 1 ) break;
                         ++mismtch2;    }
               if ( mtch2 > 0 && mismtch2 > 0 
                    && rd1[p1+mtch2-1] != rd2[p2+mtch2] )
               {    --mtch2;
                    --mismtch2;    }

               // Compute the number of matching characters if we assume a
               // double gap on read2.  Allow one extra mismatch.

               int misdmtch2 = 0, dmtch2;
               for ( dmtch2 = 0; dmtch2 < gap1 && dmtch2 + 2 < gap2; dmtch2++ )
                    if ( rd1[p1+dmtch2] != rd2[p2+dmtch2+2] )
                    {    if ( misdmtch2 == 1 ) break;
                         ++misdmtch2;    }
               if ( dmtch2 > 0 && misdmtch2 > 0 
                    && rd1[p1+dmtch2-1] != rd2[p2+dmtch2+1] )
               {    --dmtch2;
                    --misdmtch2;    }
               int mtchx = Min( mtch, Min(gap1, gap2) - 1 );
               if ( (mtch - mismtch > 0 || Min(gap1, gap2) == 1) 
                    && mtchx >= 0 && mtch - mismtch >= mtch1 - mismtch1 - 1 
                    && mtch - mismtch >= mtch2 - mismtch2 - 1
                    && mtch - mismtch >= dmtch1 - misdmtch1 - 3
                    && mtch - mismtch >= dmtch2 - misdmtch2 - 3 )
               {    if ( mtchx <= first_bug ) mismtch = 0;
                    mtch = mtchx;
                    lengths( lengths_length - 1 ) += mtch + 1;
                    p1 += mtch + 1;
                    p2 += mtch + 1;
                    gap1 -= mtch + 1;
                    gap2 -= mtch + 1;
                    errors += 1 + mismtch;    }
               else if ( mtch1 - mismtch1 > 0 
                    && mtch1 - mismtch1 >= mtch2 - mismtch2
                    && mtch1 - mismtch1 >= dmtch1 - misdmtch1 - 2
                    && mtch1 - mismtch1 >= dmtch2 - misdmtch2 - 2 )
               {    lengths(lengths_length++) = mtch1;
                    gaps(gaps_length++) = -1;
                    p1 += mtch1 + 1;
                    p2 += mtch1;
                    gap1 -= mtch1 + 1;
                    gap2 -= mtch1;
                    errors += 2 + mismtch1;    }
               else if ( mtch2 - mismtch2 > 0
                    && mtch2 - mismtch2 >= dmtch1 - misdmtch1 - 2
                    && mtch2 - mismtch2 >= dmtch2 - misdmtch2 - 2 )
               {    lengths(lengths_length++) = mtch2;
                    gaps(gaps_length++) = 1;
                    p1 += mtch2;
                    p2 += mtch2 + 1;
                    gap1 -= mtch2;
                    gap2 -= mtch2 + 1;
                    errors += 2 + mismtch2;    }
               else if ( dmtch1 - misdmtch1 > 2
                    && dmtch1 - misdmtch1 >= dmtch2 - misdmtch2 )
               {    lengths(lengths_length++) = dmtch1;
                    gaps(gaps_length++) = -2;
                    p1 += dmtch1 + 2;
                    p2 += dmtch1;
                    gap1 -= dmtch1 + 2;
                    gap2 -= dmtch1;
                    errors += 4 + misdmtch1;    }
               else if ( dmtch2 - misdmtch2 > 2 )
               {    lengths(lengths_length++) = dmtch2;
                    gaps(gaps_length++) = 2;
                    p1 += dmtch2;
                    p2 += dmtch2 + 2;
                    gap1 -= dmtch2;
                    gap2 -= dmtch2 + 2;
                    errors += 4 + misdmtch2;    }
               else if ( gap2 == 0 )
               {    p1 += gap1;
                    while ( gap1 > 19 )
                    {    gaps(gaps_length++) = -19;
                         lengths(lengths_length++) = 0;
                         gap1 -= 19;
                         errors += 38;    }
                    gaps(gaps_length++) = -gap1;
                    lengths(lengths_length++) = L;
                    errors += 2 * gap1;    
                    break;    }
               else if ( gap1 == 0 )
               {    p2 += gap2;
                    while ( gap2 > 19 )
                    {    gaps(gaps_length++) = 19;
                         lengths(lengths_length++) = 0;
                         gap2 -= 19;
                         errors += 38;    }
                    gaps(gaps_length++) = gap2;
                    lengths(lengths_length++) = L;
                    errors += 2 * gap2;    
                    break;    }
               else
               {    --gap1;
                    --gap2;
                    ++p1;
                    ++p2;
                    ++lengths(lengths_length-1);
                    ++errors;    }    }

          p1 += L;
          p2 += L;    
          if ( errors - errors_was > local_max_errors ) return;
          Assert( p1 <= rd1.size( ) ); // XXX
          Assert( p2 <= rd2.size( ) ); // XXX
          }

     // Now walk over any single base errors on the right.
     ForceAssert( p1 <= rd1.size( ) ); // XXX
     ForceAssert( p2 <= rd2.size( ) ); // XXX

     unsigned int rspace;
     int lspace;
     int errors_was = errors;

     rwalk:
     if ( errors - errors_was > local_max_errors_done ) goto done_walking;
     if ( errors - errors_was > local_max_errors ) return;
     if ( errors >= max_errors ) goto done_walking;
     rspace = Min( rd1.size( ) - p1, rd2.size( ) - p2 );
     if ( rspace >= 1 )
     {    int mtch, mtch1, mtch2;

          // Compute the number of matching characters after the first mismatch.
          // Allow one extra mismatch.
          int mismtch = 0;

          for ( mtch = 1; mtch < int(rspace); mtch++ ) 
               if ( rd1[p1+mtch] != rd2[p2+mtch] )
               {    if ( mismtch == 1 ) break;
                    ++mismtch;    }
          --mtch;

          // Compute the number of matching characters if we assume a gap on
          // rd1.  Allow one extra mismatch.
          int mismtch1 = 0;
          for ( mtch1 = 0; 
               p1 + mtch1 + 1 < rd1.size( ) && p2 + mtch1 < rd2.size( ); 
               mtch1++ )
               if ( rd1[p1+mtch1+1] != rd2[p2+mtch1] )
               {    if ( mismtch1 == 1 ) break;
                    ++mismtch1;    }
          if ( mtch1 > 0 && mismtch1 > 0 
               && rd1[p1+mtch1] != rd2[p2+mtch1-1] )
          {    --mtch1;
               --mismtch1;    }

          // Compute the number of matching characters if we assume a double 
          // gap on rd1.  Allow one extra mismatch.
          int misdmtch1 = 0, dmtch1;
          for ( dmtch1 = 0; p1 + dmtch1 + 2 < rd1.size( ) 
               && p2 + dmtch1 < rd2.size( ); dmtch1++ )
               if ( rd1[p1+dmtch1+2] != rd2[p2+dmtch1] )
               {    if ( misdmtch1 == 1 ) break;
                    ++misdmtch1;    }
          if ( dmtch1 > 0 && misdmtch1 > 0 
               && rd1[p1+dmtch1+1] != rd2[p2+dmtch1-1] )
          {    --dmtch1;
               --misdmtch1;    }

          // Compute the number of matching characters if we assume a gap on
          // read2.  Allow one extra mismatch.
          int mismtch2 = 0;
          for ( mtch2 = 0; 
               p1 + mtch2 < rd1.size( ) && p2 + mtch2 + 1 < rd2.size( ); 
               mtch2++ )
               if ( rd1[p1+mtch2] != rd2[p2+mtch2+1] )
               {    if ( mismtch2 == 1 ) break;
                    ++mismtch2;    }
          if ( mtch2 > 0 && mismtch2 > 0 
               && rd1[p1+mtch2-1] != rd2[p2+mtch2] )
          {    --mtch2;
               --mismtch2;    }

          // Compute the number of matching characters if we assume a double
          // gap on read2.  Allow one extra mismatch.
          int misdmtch2 = 0, dmtch2;
          for ( dmtch2 = 0; p1 + dmtch2 < rd1.size( ) 
               && p2 + dmtch2 + 2 < rd2.size( ); dmtch2++ )
               if ( rd1[p1+dmtch2] != rd2[p2+dmtch2+2] )
               {    if ( misdmtch2 == 1 ) break;
                    ++misdmtch2;    }
          if ( dmtch2 > 0 && misdmtch2 > 0 
               && rd1[p1+dmtch2-1] != rd2[p2+dmtch2+1] )
          {    --dmtch2;
               --misdmtch2;    }

          if ( (mtch - mismtch > 0 || rspace == 1) 
               && mtch - mismtch >= mtch1 - mismtch1 - 1 
               && mtch - mismtch >= mtch2 - mismtch2 - 1
               && mtch - mismtch >= dmtch1 - misdmtch1 - 3
               && mtch - mismtch >= dmtch2 - misdmtch2 - 3 )
          {    lengths( lengths_length - 1 ) += mtch + 1;
               p1 += mtch + 1;
               p2 += mtch + 1;
               errors += 1 + mismtch;
               goto rwalk;    }
          if ( mtch1 - mismtch1 > 0 
               && (mtch1 - mismtch1 > mtch2 - mismtch2 ||
                    (mtch1 - mismtch1 == mtch2 - mismtch2 && mismtch1 <= mismtch2))
               && mtch1 - mismtch1 >= dmtch1 - misdmtch1 - 2
               && mtch1 - mismtch1 >= dmtch2 - misdmtch2 - 2 )
          {    lengths(lengths_length++) = mtch1;
               gaps(gaps_length++) = -1;
               p1 += mtch1 + 1;
               p2 += mtch1;
               errors += 2 + mismtch1;
               goto rwalk;    }
          if ( mtch2 - mismtch2 > 0 && mtch2 - mismtch2 >= dmtch1 - misdmtch1 - 2
               && mtch2 - mismtch2 >= dmtch2 - misdmtch2 - 2 )
          {    lengths(lengths_length++) = mtch2;
               gaps(gaps_length++) = 1;
               p1 += mtch2;
               p2 += mtch2 + 1;
               errors += 2 + mismtch2;
               goto rwalk;    }
          if ( dmtch1 - misdmtch1 > 2 && dmtch1 - misdmtch1 >= dmtch2 - misdmtch2 )
          {    lengths(lengths_length++) = dmtch1;
               gaps(gaps_length++) = -2;
               p1 += dmtch1 + 2;
               p2 += dmtch1;
               errors += 4 + misdmtch1;
               goto rwalk;    }
          if ( dmtch2 - misdmtch2 > 2 )
          {    lengths(lengths_length++) = dmtch2;
               gaps(gaps_length++) = 2;
               p1 += dmtch2;
               p2 += dmtch2 + 2;
               errors += 4 + misdmtch2;
               goto rwalk;    }
          if ( int(rspace) < end_stretch*k )
          {    ++p1;
               ++p2;
               ++lengths(lengths_length-1);
               ++errors;    
               goto rwalk;    }    }

     // Now walk over any single base errors on the left.
     Assert( p1 <= rd1.size( ) ); // XXX
     Assert( p2 <= rd2.size( ) ); // XXX
     errors_was = errors;

     lwalk:
     // PRINT(errors); // XXX
     if ( errors - errors_was > local_max_errors_done ) goto done_walking;
     if ( errors - errors_was > local_max_errors ) return;
     if ( errors >= max_errors ) goto done_walking;
     lspace = Min( pos1, pos2 );
     if ( lspace >= 1 )
     {    int mtch, mtch1, dmtch1, mtch2, dmtch2;

          // Compute the number of matching characters before the first 
          // mismatch.  Allow one extra mismatch.
          int mismtch = 0;
          int mtchx = 0;
          for ( mtch = 1; mtch + 1 <= lspace; mtch++ ) 
               if ( rd1[pos1-mtch-1] != rd2[pos2-mtch-1] )
               {    if ( mismtch == 1 ) break;
                    mtchx = mtch - 1;
                    ++mismtch;    }
          --mtch;

          // Compute the number of matching characters if we assume a gap on
          // rd1.  Allow one extra mismatch.
          int mismtch1 = 0, mtch1x = 0;
          for ( mtch1 = 1; pos1 - mtch1 >= 0 && pos2-mtch1-1 >= 0; mtch1++ )
               if ( rd1[pos1-mtch1] != rd2[pos2-mtch1-1] )
               {    if ( mismtch1 == 1 ) break;
                    mtch1x = mtch1 - 1;
                    ++mismtch1;    }
          --mtch1;

          // Compute the number of matching characters if we assume a double 
          // gap on rd1.  Allow one extra mismatch.
          int misdmtch1 = 0, dmtch1x = 0;
          for ( dmtch1 = 1; pos1 - dmtch1 >= 0 && pos2-dmtch1-2 >= 0; dmtch1++ )
               if ( rd1[pos1-dmtch1] != rd2[pos2-dmtch1-2] )
               {    if ( misdmtch1 == 1 ) break;
                    dmtch1x = dmtch1 - 1;
                    ++misdmtch1;    }
          --dmtch1;

          // Compute the number of matching characters if we assume a gap on
          // read2.  Allow one extra mismatch.
          int mismtch2 = 0, mtch2x = 0;
          for ( mtch2 = 1; pos1 - mtch2 - 1 >= 0 && pos2 - mtch2 >= 0; mtch2++ )
               if ( rd1[pos1-mtch2-1] != rd2[pos2-mtch2] )
               {    if ( mismtch2 == 1 ) break;
                    mtch2x = mtch2 - 1;
                    ++mismtch2;    }
          --mtch2;

          // Compute the number of matching characters if we assume a double 
          // gap on read2.
          int misdmtch2 = 0, dmtch2x = 0;
          for ( dmtch2 = 1; pos1-dmtch2-2 >= 0 && pos2 - dmtch2 >= 0; dmtch2++ )
               if ( rd1[pos1-dmtch2-2] != rd2[pos2-dmtch2] )
               {    if ( misdmtch2 == 1 ) break;
                    dmtch2x = dmtch2 - 1;
                    ++misdmtch2;    }
          --dmtch2;

          if ( (mtch - mismtch > 0 || lspace == 1) 
               && mtch - mismtch >= mtch1 - mismtch1 - 1 
                    && mtch - mismtch >= mtch2 - mismtch2 - 1 
               && mtch - mismtch >= dmtch1 - misdmtch1 - 3 
               && mtch - mismtch >= dmtch2 - misdmtch2 - 3 )
          {    if ( mismtch == 1 && mtchx > 0 ) 
               {    mtch = mtchx;
                    mismtch = 0;    }
               lengths(lengths_start) += mtch + 1;
               pos1 -= mtch + 1;
               pos2 -= mtch + 1;
               errors += 1 + mismtch;
               goto lwalk;    }
          if ( mtch1 - mismtch1 > 0 && mtch1 - mismtch1 >= mtch2 - mismtch2 
               && mtch1 - mismtch1 >= dmtch1 - misdmtch1 - 2 &&
               mtch1 - mismtch1 >= dmtch2 - misdmtch2 - 2 )
          {    if ( mismtch1 == 1 && mtch1x > 0 ) 
               {    mtch1 = mtch1x;
                    mismtch1 = 0;    }
               lengths(--lengths_start) = mtch1;
               gaps(gaps_start) = 1;
               gaps(--gaps_start) = 0;
               pos1 -= mtch1;
               pos2 -= mtch1 + 1;
               errors += 2 + mismtch1;
               goto lwalk;    }
          if ( mtch2 - mismtch2 > 0 
               && mtch2 - mismtch2 >= dmtch1 - misdmtch1 - 2 
               && mtch2 - mismtch2 >= dmtch2 - misdmtch2 - 2 )
          {    if ( mismtch2 == 1 && mtch2x > 0 ) 
               {    mtch2 = mtch2x;
                    mismtch2 = 0;    }
               lengths(--lengths_start) = mtch2;
               gaps(gaps_start) = -1;
               gaps(--gaps_start) = 0;
               pos1 -= mtch2 + 1;
               pos2 -= mtch2;
               errors += 2 + mismtch2;
               goto lwalk;    }    
          if ( dmtch1 - misdmtch1 > 2 
               && dmtch1 - misdmtch1 >= dmtch2 - misdmtch2 )
          {    if ( misdmtch1 == 1 && dmtch1x > 0 ) 
               {    dmtch1 = dmtch1x;
                    misdmtch1 = 0;    }
               lengths(--lengths_start) = dmtch1;
               gaps(gaps_start) = 2;
               gaps(--gaps_start) = 0;
               pos1 -= dmtch1;
               pos2 -= dmtch1 + 2;
               errors += 4;
               goto lwalk;    }
          if ( dmtch2 - misdmtch2 > 2 )
          {    if ( misdmtch2 == 1 && dmtch2x > 0 ) 
               {    dmtch2 = dmtch2x;
                    misdmtch2 = 0;    }
               lengths(--lengths_start) = dmtch2;
               gaps(gaps_start) = -2;
               gaps(--gaps_start) = 0;
               pos1 -= dmtch2 + 2;
               pos2 -= dmtch2;
               errors += 4;
               goto lwalk;    }    
          if ( lspace < end_stretch*k )
          {    --pos1;
               --pos2;
               if ( gaps(gaps_start) == 0 ) ++lengths(lengths_start);
               else
               {    gaps(--gaps_start) = 0;
                    lengths(--lengths_start) = 1;    }
               ++errors;
               goto lwalk;    }    }

     done_walking:
     Assert( gaps_length - gaps_start == lengths_length - lengths_start );

     Assert( p1 <= rd1.size( ) ); // XXX
     Assert( p2 <= rd2.size( ) ); // XXX

     Assert( gaps_start > 1 ) ; 
     Assert( gaps_length < 2 * halfgaps - 1 ) ; 

     avector<int> gaps2(2 * halfgaps), lengths2(2 * halfgaps);
     int l = 0;
     for ( int j = gaps_start; j < gaps_length; j++ )
     {    if ( j == gaps_start || gaps(j) != 0 )
          {    gaps2(l) = gaps(j);
               lengths2(l) = lengths(j);
               ++l;    }
          else lengths2( l - 1 ) += lengths(j);    }

     Setpos1(pos1);
     Setpos2(pos2);
     SetNblocks(l);
     for ( int i = 0; i < l; i++ )
     {    SetGap( i, gaps2(i) );
          SetLength( i, lengths2(i) );    }
     errors_found = errors;    }
