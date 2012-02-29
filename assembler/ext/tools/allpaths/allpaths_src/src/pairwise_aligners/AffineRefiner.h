// Copyright (c) 2004  The Broad Institute/Massachusetts Institute of Technology

#ifndef AFFINE_REFINER_H
#define AFFINE_REFINER_H

#include "Alignment.h"
#include "Basevector.h"
#include "Vec.h"

// AffineRefiner

// This class takes an alignment, finds its gappy regions, and
// performs an affine SW alignment on those regions, pasting the
// results back into the original alignment.

class AffineRefiner
{
 public:
  // A length in the alignment must be at least this many bases long
  // to be considered an anchor.  Only gappy regions *between* anchors
  // will be patched.
  AffineRefiner( const int minAnchorLength = 12 )
    : m_minAnchorLength( minAnchorLength )
  {}

  // Refine the given alignment between the given sequences.
  void
  RefineAlign( align &theAlign,
               const basevector &bases1, 
               const basevector &bases2 );
  
 private:
  // Check whether the given block of the alignment begins or ends a gappy region.
  void CheckBlock( const int currentBlock, 
                   const align &theAlign,
                   const basevector &bases1, const basevector &bases2 );

  bool BlockIsAnchor( const int currentBlock,
                      const align &theAlign,
                      const basevector &bases1, const basevector &bases2 );
  
  // Patch the given region with a SW affine alignment.
  void PatchAlignment( const align &theAlign, 
                       const int startAnchor, 
                       const int stopAnchor, 
                       const basevector &bases1, const basevector &bases2 );

  int m_minAnchorLength;

  vec<int> m_newAlignGaps;
  vec<int> m_newAlignLens;
  align    m_newAlign;

  int m_posOn1, m_posOn2;
  int m_lastAnchor;
  int m_sumGaps;
  int m_sumLengths;
  bool m_foundAnchor;

  basevector m_chunk1, m_chunk2;
  alignment  m_patchAlignment;
  align      m_patchAlign;
};



#endif
