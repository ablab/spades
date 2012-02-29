// Copyright (c) 2004  The Broad Institute/Massachusetts Institute of Technology

#include "pairwise_aligners/AffineRefiner.h"

#include "pairwise_aligners/SmithWatAffine.h"

void
AffineRefiner::RefineAlign( align &theAlign, 
                            const basevector &bases1, 
                            const basevector &bases2 )
{
  /* FOR DEBUGGING
  cout << "old alignment:" << endl;
  {
    int pos1 = theAlign.pos1();
    int pos2 = theAlign.pos2();
    for ( int blockNum = 0; blockNum < theAlign.Nblocks(); ++blockNum )
    {
      int gap = theAlign.Gaps( blockNum );
      int len = theAlign.Lengths( blockNum );
      if ( gap < 0 )
        pos1 -= gap;
      else
        pos2 += gap;
      PRINT4( gap, pos1, pos2, len );
      pos1 += len;
      pos2 += len;
    }
  }
  PrintVisualAlignment( True, cout, bases1, bases2, theAlign );
  */


  // Find the gaps and lengths for the new align.
  m_newAlignGaps.clear();
  m_newAlignLens.clear();

  m_sumGaps = 0;
  m_sumLengths = 0;
  m_lastAnchor = -1;
  m_posOn1 = theAlign.pos1();
  m_posOn2 = theAlign.pos2();
  
  for ( int blockNum = 0; blockNum < theAlign.Nblocks(); ++blockNum )
    this->CheckBlock( blockNum, theAlign,
                      bases1, bases2 );
  
  for ( int copyBlockNum = m_lastAnchor + 1;
        copyBlockNum < theAlign.Nblocks(); ++copyBlockNum )
  {
    m_newAlignGaps.push_back( theAlign.Gaps( copyBlockNum ) );
    m_newAlignLens.push_back( theAlign.Lengths( copyBlockNum ) );
  }

  // Build the new align.
  ForceAssertEq( m_newAlignGaps.size(), m_newAlignLens.size() );
  m_newAlign.Setpos1( theAlign.pos1() );
  m_newAlign.Setpos2( theAlign.pos2() );
  m_newAlign.SetNblocks( m_newAlignGaps.size() );
  
  if ( m_newAlign.Nblocks() > 0 )
  {
    int firstGap = m_newAlignGaps.front();
    if ( firstGap < 0 )
      m_newAlign.Setpos1( m_newAlign.pos1() - firstGap );
    if ( firstGap > 0 )
      m_newAlign.Setpos2( m_newAlign.pos2() + firstGap );
    m_newAlignGaps.front() = 0;
  }

  for ( int ii = 0; ii < m_newAlign.Nblocks(); ++ii )
  {
    m_newAlign.SetGap( ii, m_newAlignGaps[ii] );
    m_newAlign.SetLength( ii, m_newAlignLens[ii] );
  }

  // Copy it to the input align.
  theAlign = m_newAlign;

  /* FOR DEBUGGING
  cout << "new alignment:" << endl;
  {
    int pos1 = theAlign.pos1();
    int pos2 = theAlign.pos2();
    for ( int blockNum = 0; blockNum < theAlign.Nblocks(); ++blockNum )
    {
      int gap = theAlign.Gaps( blockNum );
      int len = theAlign.Lengths( blockNum );
      if ( gap < 0 )
        pos1 -= gap;
      else
        pos2 += gap;
      PRINT4( gap, pos1, pos2, len );
      pos1 += len;
      pos2 += len;
    }
  }
  PrintVisualAlignment( True, cout, bases1, bases2, theAlign );
  */
}


bool
AffineRefiner::BlockIsAnchor( const int currentBlock,
                              const align &theAlign,
                              const basevector &bases1,
                              const basevector &bases2 )
{
  int len = theAlign.Lengths( currentBlock );
  if ( len < m_minAnchorLength ) 
    return false;

  int perfectRunLength = 0;
  for ( int ii = 0; ii < len; ++ii )
    if ( bases1[ m_posOn1 + ii ] != bases2[ m_posOn2 + ii ] )
      perfectRunLength = 0;
    else if ( ++perfectRunLength >= m_minAnchorLength )
      return true;

  return false;
}

    
void 
AffineRefiner::CheckBlock( const int currentBlock, 
                           const align &theAlign, 
                           const basevector &bases1, const basevector &bases2 )
{
  int gap = theAlign.Gaps( currentBlock );
  if ( gap < 0 )
  {
    m_posOn1 -= gap;
    m_sumGaps -= gap;
  }
  else
  {
    m_posOn2 += gap;
    m_sumGaps += gap;
  }

  int mismatches, maxPerfectMatch;
  if ( this->BlockIsAnchor( currentBlock, theAlign, bases1, bases2 ) )
  {
    // If the amount of gaps between this anchor and the last is at
    // least half the amount of lengths AND there's more than one gap
    // between this anchor and the last, patch it.
    if ( m_sumGaps >= m_sumLengths/2  &&
         m_lastAnchor >= 0 &&
         currentBlock - m_lastAnchor > 1 )
      this->PatchAlignment( theAlign,
                            m_lastAnchor, currentBlock,
                            bases1, bases2 );

    // Otherwise, just copy the intervening blocks from the old align.
    else
    {
      for ( int copyBlockNum = m_lastAnchor + 1;
            copyBlockNum <= currentBlock; ++copyBlockNum )
      {
        m_newAlignGaps.push_back( theAlign.Gaps( copyBlockNum ) );
        m_newAlignLens.push_back( theAlign.Lengths( copyBlockNum ) );
      }
    }

    m_lastAnchor = currentBlock;
    m_sumGaps = 0;
    m_sumLengths = 0;
  }

  // If this block isn't an anchor, proceed.
  else
    m_sumLengths += theAlign.Lengths( currentBlock );
  
  m_posOn1 += theAlign.Lengths( currentBlock );
  m_posOn2 += theAlign.Lengths( currentBlock );
}


void
AffineRefiner::PatchAlignment( const align &theAlign, 
                               const int startAnchor, 
                               const int stopAnchor, 
                               const basevector &bases1, const basevector &bases2 )
{
  // Find where the last block added to the new alignment ends, which
  // is where our patch begins.
  int start1 = theAlign.pos1();
  int start2 = theAlign.pos2();

  for ( int newBlockNum = 0; newBlockNum < (int) m_newAlignGaps.size(); ++newBlockNum )
  {
    int gap = m_newAlignGaps[ newBlockNum ];
    if ( gap < 0 )
      start1 -= gap;
    else
      start2 += gap;

    int len = m_newAlignLens[ newBlockNum ];
    start1 += len;
    start2 += len;
  }

  // Find the length of the gappy region we're patching.
  int len1 = 0;
  int len2 = 0;
  
  for ( int blockNum = startAnchor + 1; blockNum <= stopAnchor; ++blockNum )
  {
    int gap = theAlign.Gaps( blockNum );
    if ( gap < 0)
      len1 -= gap;
    else
      len2 += gap;

    if ( blockNum < stopAnchor )
    {
      int len = theAlign.Lengths( blockNum );
      len1 += len;
      len2 += len;
    }
  }

  m_chunk1.SetToSubOf( bases1, start1, len1 );
  m_chunk2.SetToSubOf( bases2, start2, len2 );

  int bestOffset;
  alignment m_patchAlignment;
  const bool penalizeLeftGap = true;
  const bool penalizeRightGap = true;
  const int mismatchPenalty = 3;
  const int gapOpenPenalty = 6;
  const int gapExtendPenalty = 1;
  int score = SmithWatAffine( m_chunk1, m_chunk2, m_patchAlignment,
                              penalizeLeftGap, penalizeRightGap,
                              mismatchPenalty, gapOpenPenalty, gapExtendPenalty );
  align m_patchAlign( m_patchAlignment );

  /* FOR DEBUGGING
  cout << "chunk alignment:" << endl;
  PRINT2( start1, len1 );
  m_chunk1.Print( cout, "chunk1" );
  PRINT2( start2, len2 );
  m_chunk2.Print( cout, "chunk2" );
  PrintVisualAlignment( True, cout, m_chunk1, m_chunk2, m_patchAlign );
  */

  ForceAssertGt( m_patchAlign.Nblocks(), 0 );

  // Add a gap if the sequences don't line up on the left.
  ForceAssert( ! ( m_patchAlign.pos1() > 0 && m_patchAlign.pos2() > 0 ) );

  if ( m_patchAlign.pos1() > 0 )
  {
    m_newAlignGaps.push_back( -m_patchAlign.pos1() );
    m_newAlignLens.push_back( 0 );
  }
  if ( m_patchAlign.pos2() > 0 )
  {
    m_newAlignGaps.push_back( m_patchAlign.pos2() );
    m_newAlignLens.push_back( 0 );
  }

  // Add the patch alignment.
  int blockIdx = 0;

  // If the first gap is zero (which it should be), add the first
  // length of the patch to the last length of the new alignment and
  // move on to the next block.
  if ( m_patchAlign.Gaps( 0 ) == 0 )
  {
    m_newAlignLens.back() += m_patchAlign.Lengths( 0 );
    ++blockIdx;
  }

  for ( ; blockIdx < m_patchAlign.Nblocks(); ++blockIdx )
  {
    m_newAlignGaps.push_back( m_patchAlign.Gaps( blockIdx ) );
    m_newAlignLens.push_back( m_patchAlign.Lengths( blockIdx ) );
  }

  // Add a trailing gap if the sequences don't line up on the right.
  int tail1 = (int)m_chunk1.size() - m_patchAlign.Pos1();
  int tail2 = (int)m_chunk2.size() - m_patchAlign.Pos2();

  ForceAssert( ! ( tail1 > 0 && tail2 > 0 ) );

  if ( tail1 > 0 )
  {
    m_newAlignGaps.push_back( -tail1 );
    m_newAlignLens.push_back( 0 );
  }
  if ( tail2 > 0 )
  {
    m_newAlignGaps.push_back( tail2 );
    m_newAlignLens.push_back( 0 );
  }

  // Add the anchoring length segment.
  m_newAlignLens.back() += theAlign.Lengths( stopAnchor );
}

