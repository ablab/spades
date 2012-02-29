///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef FEUDAL_BITVEC_H_
#define FEUDAL_BITVEC_H_

#include "feudal/FieldVec.h"
#include "feudal/MasterVec.h"
#include "feudal/Mempool.h"
#include "feudal/BinaryStreamTraits.h"
#include "system/Assert.h"
#include <algorithm>
#include <cstring>
#include <ostream>

/// A feudal vector that has bits as values.
class BitVec : public FieldVec<1, MempoolAllocator<unsigned char> >
{
public:
    typedef FieldVec<1, MempoolAllocator<unsigned char> > Base;
    typedef MempoolAllocator<unsigned char> Alloc;

    // Constructors
    BitVec() : Base() {}
    BitVec(Alloc const& alloc) : Base(alloc) {}
    BitVec(Base::size_type n, Base::value_type exemplar = 0,
            Base::size_type cap = 0, Alloc const& alloc = Alloc())
    : Base(n, exemplar, cap, alloc) {}

    BitVec& operator|=( BitVec const& bv );
    BitVec& operator&=( BitVec const& bv );
    BitVec& operator^=( BitVec const& bv );

    // Wrappers
    void Set(size_type i, value_type bit) { set(i, bit); }

    BitVec& ReverseMe() { std::reverse(begin(),end()); return *this; }

    /// Return the index of the next element after i that is different from
    /// the ith element.  Returns size() if there is no such element.
    size_type NextDiff(Base::size_type i) const
    { Base::value_type val = (*this)[i];
      for ( ++i; i < size(); ++i )
        if ( (*this)[i] != val) break;
      return i; }

    BitVec& Zero() { size_type n = size(); clear().resize(n); return *this; }

    BitVec& invert();

    // Set *this to the length len sub-bitvector of src, starting at position
    // start_pos.  The case where this == &src is allowed.
    BitVec& SetToSubOf(const BitVec& src, size_type start_pos, size_type len);


    // Set bits x (start <= x < stop) to the same value "bit".
    BitVec& Set( size_type start, size_type stop, value_type bit )
    { AssertLe(stop,size()); AssertLe(start,stop);
      iterator end(begin(stop));
      for ( iterator itr(begin(start)); itr != end; ++itr )
        itr.set(bit);
      return *this; }

    // Return a sum of the bits
    size_type Sum() const
    { size_type sum = 0;
      for ( const_iterator itr(begin()), stop(end()); itr != stop; ++itr )
        sum += *itr;
      return sum; }

    void PrintFastaStyle(ostream& out, const String& id) const;

    friend ostream& operator<<( ostream& s, const BitVec& bv )
    { BitVec::const_iterator stop(bv.end());
      for ( BitVec::const_iterator itr(bv.begin()); itr != stop; ++itr )
        s << (*itr ? '1' : '0');
      return s; }
};

SELF_SERIALIZABLE(BitVec);

namespace std
{
template<> inline void iter_swap( BitVec::iterator itr1,
                                  BitVec::iterator itr2 )
{
    BitVec::value_type tmp = *itr1;
    itr1.set(*itr2);
    itr2.set(tmp);
}

template<> inline void iter_swap( BitVec::reverse_iterator itr1,
                                  BitVec::reverse_iterator itr2 )
{
    BitVec::value_type tmp = *itr1;
    itr1.set(*itr2);
    itr2.set(tmp);
}
}

inline void swap( BitVec& bv1, BitVec& bv2 )
{
    bv1.swap(bv2);
}

typedef MasterVec<BitVec> vecbitvector;

float Coverage(const vecbitvector& v);

vecbitvector& operator |=( vecbitvector& v1, const vecbitvector& v2 );
vecbitvector& operator &=( vecbitvector& v1, const vecbitvector& v2 );
vecbitvector& operator ^=( vecbitvector& v1, const vecbitvector& v2 );

#endif /* FEUDAL_BITVEC_H_ */
