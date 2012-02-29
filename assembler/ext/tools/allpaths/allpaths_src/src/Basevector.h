///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BASEVECTOR_H_
#define BASEVECTOR_H_

#include "feudal/BaseVec.h"
#include "feudal/MasterVec.h"
#include <cstddef>

typedef BaseVec basevector;
typedef BaseVec bvec;

typedef MasterVec<BaseVec> BaseVecVec;
typedef BaseVecVec vecbasevector;
typedef BaseVecVec vecbvec;

typedef OuterVec< OuterVec<bvec,MempoolAllocator<unsigned char> >,
                  MempoolOwner<unsigned char> > bvec3;

typedef bvec kmer_t;
typedef bvec read_t;
typedef bvec genome_part_t;
typedef vecbvec reads_t;
typedef vecbvec genome_t;

void ReverseComplement( vecbvec& s );

/// A class that iterates through a vecbvec as if it were concatenated into
/// a single bvec.
class FlatIterator
: public std::iterator<std::input_iterator_tag,unsigned char,
                       std::ptrdiff_t, void, unsigned char>
{
public:
    FlatIterator( vecbvec& vbv, unsigned long posn = 0 )
    : mVBV(vbv), mOuterPosn(posn), mInnerPosn(INIT_IPOSN)
    { advance(); }

    // compiler-supplied copying and destructor are OK

    FlatIterator& operator++() { advance(); return *this; }

    FlatIterator operator++( int )
    { FlatIterator tmp(*this); advance(); return tmp; }

    unsigned char operator*() const { return mVBV[mOuterPosn][mInnerPosn]; }

    friend bool operator==( FlatIterator const& fi1, FlatIterator const& fi2 )
    { return &fi1.mVBV == &fi2.mVBV &&
             fi1.mOuterPosn == fi2.mOuterPosn &&
             fi1.mInnerPosn == fi2.mInnerPosn; }

    friend bool operator!=( FlatIterator const& fi1, FlatIterator const& fi2 )
    { return !(fi1==fi2); }

private:
    void advance()
    { while ( mOuterPosn < mVBV.size() &&
              ++mInnerPosn >= mVBV[mOuterPosn].size() )
      { ++mOuterPosn; mInnerPosn = INIT_IPOSN; } }

    vecbvec& mVBV;
    vecbvec::size_type mOuterPosn;
    bvec::size_type mInnerPosn;

    static bvec::size_type const INIT_IPOSN = static_cast<bvec::size_type>(0)-1;
};

#endif // BASEVECTOR_H_
