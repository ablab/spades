///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef COMPRESSEDSEQUENCE_H_
#define COMPRESSEDSEQUENCE_H_

#include "Vec.h"
#include "Basevector.h"
#include "Bitvector.h"
#include "feudal/FieldVec.h"
#include "feudal/MasterVec.h"
#include "feudal/Mempool.h"
#include <algorithm>
#include <cstring>

// This class stores constant sequences of [ACTGN] in some fraction of the space
// required for full text, with reasonably speedy expansion and compaction.
// For simplicity's sake, they are more or less required to be const.

class CompressedSequence : public FieldVec<4, MempoolAllocator<unsigned char> >
{
public:
    typedef allocator_type Alloc;
    typedef FieldVec<4, MempoolAllocator<unsigned char> > BaseT;

    CompressedSequence() {}
    CompressedSequence( Alloc const& alloc ) : BaseT(alloc) {}

    // Copy constructor
    CompressedSequence( CompressedSequence const& cs ) : BaseT(cs) {}

    CompressedSequence( char const* str ) { assignChars(str,str+strlen(str)); }

    CompressedSequence( char const* start, char const* end )
    { assignChars(start,end); }

    CompressedSequence( const vec<char>& cv )
    { char const* buf = &*cv.begin();
      char const* end = buf+cv.size();
      assignChars(buf,end); }

    CompressedSequence( const basevector& bv )
    : BaseT(bv.begin(),bv.end(),Base::val2Bits) {}

    // compiler-supplied copy-assignment and destructor are OK

    void ReverseComplement();

    vec<char> asVecChar() const
    { vec<char> result; asVecChar(result); return result; }

    vec<char> SubAsVecChar( int begin, int end ) const
    { vec<char> result; SubAsVecChar(result,begin,end); return result; }

    basevector asBasevector() const
    { bvec result; asBasevector(result); return result; }

    String asString() const
    { String result(size(),'X');
      std::transform(begin(),end(),result.begin(),GeneralizedBase::bits2Char);
      return result; }

    // pass-by-reference versions of "as" methods
    void asVecChar( vec<char> &cv ) const
    { cv.clear(); cv.resize(size());
      std::transform(begin(),end(),cv.begin(),GeneralizedBase::bits2Char); }

    void SubAsVecChar( vec<char> &cv, int start, int stop ) const
    { cv.clear(); cv.resize(stop-start);
      std::transform(begin(start),begin(stop),cv.begin(),GeneralizedBase::bits2Char); }

    void asBasevector( basevector &bv ) const;
    void getAmbBases( bitvector &bitv ) const;

    void assignChars( char const* begin, char const* end );
};

SELF_SERIALIZABLE(CompressedSequence);
typedef MasterVec<CompressedSequence> veccompseq;

#endif
