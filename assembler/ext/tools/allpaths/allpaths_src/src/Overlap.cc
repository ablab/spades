// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#include "Alignment.h"
#include "Basevector.h"
#include "PackAlign.h"
#include "ShortVector.h"
#include "Overlap.h"

// Note: There are two versions of EstimatedOverlap, one for an alignment, and one
// for an align.  Probably the goal should be to get rid of the one for an
// alignment.

// ================================================================================
//
// EstimatedOverlap: return the estimated overlap between two reads, or zero
// if it is unclear that they overlap at all.
//
// Note the change from pos1 == 0 to pos1 <= 1, and likewise for pos2.  I don't
// know why these were needed, but otherwise some excellent alignments are missed.
// The same change occurs in AlignFromMutmers.cc.
//
// Only the lengths of rd1 and rd2 are used.
//
// ================================================================================

int EstimatedOverlap( const alignment& a,
     const basevector& rd1, const basevector& rd2 )
{    avector<int> gaps, lengths;
     int pos1, pos2, errors;
     a.Unpack( pos1, pos2, errors, gaps, lengths );
     int p1 = pos1, p2 = pos2, total_len = 0;
     if ( gaps.length != 0 && gaps(0) != 0 )
          cout << "In EstimatedOverlap: weird alignment observed, "
               << "starting with a gap of " << gaps(0) << ".\n";
     for ( unsigned int j = 0; j < gaps.length; j++ )
     {    if ( gaps(j) > 0 ) p2 += gaps(j);
          else if ( gaps(j) < 0 ) p1 -= gaps(j);
          p1 += lengths(j); p2 += lengths(j);
          total_len += lengths(j);    }
     // The following used to have pos1 == 0, pos2 == 0.
     // Also it used to require that errors <= total_len / 10.
     if ( ( (pos1 <= 1 && (unsigned int) p2 == rd2.size( )) 
          || (pos2 <= 1 && (unsigned int) p1 == rd1.size( )) 
          || (pos1 <= 1 && (unsigned int) p1 == rd1.size( ))
          || (pos2 <= 1 && (unsigned int) p2 == rd2.size( )) ) )
          return total_len;
     else return 0;    }

int EstimatedOverlap( const align& a, const basevector& rd1, const basevector& rd2 )
{    const avector<int>& gaps = a.Gaps( ); 
     const avector<int>& lengths = a.Lengths( );
     int pos1 = a.pos1( );
     int pos2 = a.pos2( );
     int nblocks = a.Nblocks( );
     int p1 = pos1, p2 = pos2, total_len = 0;
     if ( nblocks != 0 && gaps(0) != 0 )
          cout << "In EstimatedOverlap: weird alignment observed, "
               << "starting with a gap of " << gaps(0) << ".\n";
     for ( int j = 0; j < nblocks; j++ )
     {    if ( gaps(j) > 0 ) p2 += gaps(j);
          else if ( gaps(j) < 0 ) p1 -= gaps(j);
          p1 += lengths(j); p2 += lengths(j);
          total_len += lengths(j);    }
     if ( ( (pos1 <= 1 && (unsigned int) p2 == rd2.size( )) 
          || (pos2 <= 1 && (unsigned int) p1 == rd1.size( )) 
          || (pos1 <= 1 && (unsigned int) p1 == rd1.size( ))
          || (pos2 <= 1 && (unsigned int) p2 == rd2.size( )) ) )
          return total_len;
     else return 0;    }
