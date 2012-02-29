///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CompressedSequence.h"
#include "feudal/FieldVecDefs.h"
#include "feudal/OuterVecDefs.h"

template class FieldVec< 4, MempoolAllocator<unsigned char> >;
template class OuterVec<CompressedSequence>;

void CompressedSequence::ReverseComplement()
{
    iterator head = begin();
    iterator tail = end();
    while (head != tail)
    {
        value_type tmp = GeneralizedBase::fromBits(*--tail).complement().bits();
        if (head == tail)
        {
            head.set(tmp);
            break;
        }
        tail.set(GeneralizedBase::fromBits(*head).complement().bits());
        head.set(tmp);
        ++head;
    }
}

void CompressedSequence::asBasevector( basevector &bv ) const
{
    bv.clear().reserve(size());
    for ( const_iterator itr(begin()), stop(end()); itr != stop; ++itr )
        bv.push_back(GeneralizedBase::bits2Val(*itr));
}

void CompressedSequence::getAmbBases( bitvector &bitv ) const
{
    bitv.clear().reserve(size());
    for ( const_iterator itr(begin()), stop(end()); itr != stop; ++itr )
        bitv.push_back(GeneralizedBase::bits2Ambig(*itr));
}

void CompressedSequence::assignChars( char const* begin, char const* end )
{
    clear().reserve(end-begin);
    while ( begin != end )
    {
        char chr = *begin++;
        if ( chr != '*' )
        {
            if ( GeneralizedBase::isGeneralizedBase(chr) )
                push_back(GeneralizedBase::char2Bits(chr));
            else
                push_back(GeneralizedBase::N.bits());
        }
    }
}
