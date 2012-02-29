///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file KMer.h
 * \author tsharpe
 * \date Sep 30, 2011
 *
 * \brief
 */
#ifndef KMER_H_
#define KMER_H_

#include "Compare.h"
#include "dna/Bases.h"
#include "feudal/BinaryStream.h"
#include "feudal/Iterator.h"
#include "system/StaticAssert.h"
#include <iterator>
#include <limits>
#include <ostream>

enum CanonicalForm { FWD, REV, PALINDROME };

inline std::ostream& operator<<( std::ostream& os, CanonicalForm form )
{ return os << "+-|"[form]; }


// S must be some type of unsigned integer
template <unsigned K, class S = unsigned long>
class KMer
{
public:
    typedef S storage_type;
    typedef unsigned char value_type;
    typedef unsigned size_type;
    typedef int difference_type;
    typedef std::iterator<std::random_access_iterator_tag,
                          value_type,
                          difference_type,
                          void,
                          value_type> ItrTag;

    struct NopMapper
    { value_type operator()( value_type val ) const
      { AssertLt(val,4u); return val; } };

    KMer()
    { storage_type* end(mVal + STORAGE_UNITS_PER_KMER);
      for ( storage_type* itr = mVal; itr != end; ++itr )
        *itr = 0; }


    // *Itr is a base code
    template <class Itr>
    explicit KMer( Itr start, typename Itr::iterator_category =
                                  typename Itr::iterator_category() )
    { assign(start,NopMapper()); }

    // Mapper(*Itr) is a base code
    template <class Itr, class Mapper>
    KMer( Itr start, Mapper mapper, typename Itr::iterator_category =
                                        typename Itr::iterator_category() )
    { assign(start,mapper); }

    KMer( char const* start )
    { assign(start,CharToBaseMapper()); }

    // compiler-supplied copying and destructor is OK

    class iterator
    : public ItrTag,
      public IteratorBase<iterator,size_type,difference_type>
    {
    public:
        iterator() : mpKMer(0) {}
        iterator( KMer const* pKMer, size_type pos )
        : IteratorBase<iterator,size_type,difference_type>(pos),
          mpKMer(pKMer) {}

        // compiler-supplied copying and destructor are OK

        value_type operator*() const { return mpKMer[0][this->pos()]; }

        value_type operator[]( difference_type idx ) const
        { return mpKMer[0][this->pos()+idx]; }

    protected:
        KMer const* mpKMer;
    };

    class rc_iterator
    : public ItrTag,
      public IteratorBase<rc_iterator,size_type,difference_type>
    {
    public:
        rc_iterator() : mpKMer(0), mLast(~0u) {}
        rc_iterator( KMer const* pKMer, size_type pos )
        : IteratorBase<rc_iterator,size_type,difference_type>(pos),
          mpKMer(pKMer), mLast(pKMer->size()-1u) {}

        // compiler-supplied copying and destructor are OK

        value_type operator*() const
        { return GetComplementaryBase(mpKMer[0][mLast-this->pos()]); }

        value_type operator[]( difference_type idx ) const
        { return GetComplementaryBase(mpKMer[0][mLast-(this->pos()+idx)]); }

    protected:
        KMer const* mpKMer;
        size_type mLast;
    };

    iterator begin() const { return iterator(this,0u); }
    iterator end() const { return iterator(this,size()); }
    rc_iterator rcbegin() const { return rc_iterator(this,0u); }
    rc_iterator rcend() const { return rc_iterator(this,size()); }

    static size_type getK() { return K; }
    size_type size() const { return K; }

    value_type at( size_type idx ) const
    { ForceAssertLt(idx,size()); return operator[](idx); }

    value_type operator[]( size_type idx ) const
    { AssertLt(idx,size()); return (data(idx) >> shift(idx)) & BASE_MASK; }

    value_type front() const { return operator[](0u); }
    value_type back() const { return operator[](K-1u); }

    KMer& set( size_type idx, value_type val )
    { AssertLt(idx,size()); AssertLt(val,4u);
      storage_type& datum = data(idx);
      size_type shiftCount = shift(idx);
      datum ^= (datum ^ (storage_type(val&BASE_MASK) << shiftCount))
                      & (storage_type(BASE_MASK) << shiftCount);
      return *this; }

    KMer& setFront( value_type val ) { set(0u,val); return *this; }
    KMer& setBack( value_type val ) { set(K-1u,val); return *this; }

    // *Itr is a base code
    template <class Itr>
    KMer& assign( Itr const& itr )
    { return assign(itr,NopMapper()); }

    // Mapper(*Itr) is a base code
    template <class Itr, class Mapper>
    KMer& assign( Itr itr, Mapper mapper )
    { storage_type* oItr(mVal); storage_type val = 0;
      for ( size_type idx = 1u; idx <= K; ++idx )
      { val = (val << BITS_PER_BASE) | mapper(*itr); ++itr;
        if ( !(idx%BASES_PER_STORAGE_UNIT) ) *oItr++ = val; }
      if ( UNUSED_TRAILING_BITS ) *oItr = val << UNUSED_TRAILING_BITS;
      return *this; }

    CanonicalForm getCanonicalForm() const
    { iterator fwd(begin()), rev(end());
      while ( fwd <= --rev )
      { value_type f = *fwd;
        value_type r = GetComplementaryBase(*rev);
        if ( f < r ) return FWD;
        if ( r < f ) return REV;
        ++fwd; }
      return PALINDROME; }

    bool isPalindrome() const
    { return !(K&1) && getCanonicalForm() == PALINDROME; }

    KMer& toPredecessor( value_type val )
    { AssertLt(val,4u);
      storage_type prev(storage_type(val) << SHIFT_LSBASE_TO_MSBASE);
      if ( STORAGE_UNITS_PER_KMER == 1 )
        mVal[0] = ((mVal[0] >> BITS_PER_BASE) & ~FINAL_MASK) | prev;
      else
      { storage_type* end(mVal+STORAGE_UNITS_PER_KMER);
        for ( storage_type* itr(mVal); itr != end; ++itr )
        { storage_type cur = *itr;
          *itr = (cur >> BITS_PER_BASE) | prev;
          prev = cur << SHIFT_LSBASE_TO_MSBASE; }
        end[-1] &= ~FINAL_MASK; }
      return *this; }

    KMer& toSuccessor( value_type val )
    { AssertLt(val,4u);
      storage_type prev((storage_type(val)&BASE_MASK) << UNUSED_TRAILING_BITS);
      if ( STORAGE_UNITS_PER_KMER == 1 )
        mVal[0] = (mVal[0] << BITS_PER_BASE) | prev;
      else
      { storage_type* beg(mVal);
        storage_type* itr(mVal+STORAGE_UNITS_PER_KMER);
        while ( itr-- != beg )
        { storage_type cur = *itr;
          *itr = (cur << BITS_PER_BASE) | prev;
          prev = cur >> SHIFT_LSBASE_TO_MSBASE; } }
      return *this; }

    KMer& rc()
    { if ( !UNUSED_TRAILING_BITS && sizeof(mVal)%2 == 0 )
      { unsigned char* itr1 = reinterpret_cast<unsigned char*>(mVal);
        unsigned char* itr2 = itr1+sizeof(mVal);
        while ( itr1 < itr2 )
        { unsigned char tmp = Base::rcByte(*--itr2);
          *itr2 = Base::rcByte(*itr1);
          *itr1++ = tmp; } }
      else
      { KMer that(*this);
        storage_type* oItr(mVal);
        storage_type val = 0;
        storage_type* iItr(that.mVal+STORAGE_UNITS_PER_KMER-1);
        storage_type val2 = *iItr >> UNUSED_TRAILING_BITS;
        size_type idx2 = UNUSED_TRAILING_BITS/BITS_PER_BASE;
        for ( size_type idx = 1u; idx <= K; ++idx )
        { val = (val << BITS_PER_BASE) | GetComplementaryBase(val2 & BASE_MASK);
          if ( !(idx%BASES_PER_STORAGE_UNIT) ) *oItr++ = val;
          val2 >>= 2;
          if ( !(++idx2%BASES_PER_STORAGE_UNIT) ) val2 = *--iItr; }
        if ( UNUSED_TRAILING_BITS ) *oItr = val << UNUSED_TRAILING_BITS; }
      return *this; }

    unsigned long hash() const
    { unsigned long result = 14695981039346656037ul;
      unsigned char const* itr(reinterpret_cast<unsigned char const*>(mVal));
      unsigned char const* end(itr+sizeof(mVal));
      for ( ; itr != end; ++itr )
          result = 1099511628211ul*(result ^ *itr); // FNV-1a algorithm
      return result; }

    struct Hasher : public std::unary_function<KMer,unsigned long>
    { unsigned long operator()( KMer const& kmer ) const { return kmer.hash(); } };

    template <class Itr, class OItr>
    static void kmerize( Itr beg, Itr const& end, OItr out )
    { if ( std::distance(beg,end) >= K )
      { KMer kkk(beg);
        Itr itr(beg+K);
        *out = kkk.getCanonicalForm()==FWD ? kkk : KMer(kkk).rc();
        while ( itr != end )
        { kkk.toSuccessor(*itr); ++itr; ++out;
          *out = kkk.getCanonicalForm()==FWD ? kkk : KMer(kkk).rc(); }
        ++out; } }

    template <class Itr, class OItr>
    static void kmerizeNonCanonically( Itr beg, Itr const& end, OItr out )
    { if ( std::distance(beg,end) >= K )
      { KMer kkk(beg);
        Itr itr(beg+K);
        *out = kkk;
        while ( itr != end )
        { kkk.toSuccessor(*itr); ++itr; ++out; *out = kkk; }
        ++out; } }

    friend bool operator<( KMer const& k1, KMer const& k2 )
    { return compare(k1,k2) < 0; }

    friend bool operator<=( KMer const& k1, KMer const& k2 )
    { return compare(k1,k2) <= 0; }

    friend bool operator>( KMer const& k1, KMer const& k2 )
    { return compare(k1,k2) > 0; }

    friend bool operator>=( KMer const& k1, KMer const& k2 )
    { return compare(k1,k2) >= 0; }

    friend bool operator==( KMer const& k1, KMer const& k2 )
    { return compare(k1,k2) == 0; }

    friend bool operator!=( KMer const& k1, KMer const& k2 )
    { return compare(k1,k2) != 0; }

    friend int compare( KMer const& k1, KMer const& k2 )
    { int result;
      storage_type const* end(k1.mVal + STORAGE_UNITS_PER_KMER);
      storage_type const* itr2(k2.mVal);
      for ( storage_type const* itr1 = k1.mVal; itr1 != end; ++itr1, ++itr2 )
        if ( (result = compare(*itr1,*itr2)) ) break;
      return result; }

    friend std::ostream& operator<<( std::ostream& os, KMer const& kmer )
    { for ( iterator itr(kmer.begin()), end(kmer.end()); itr != end; ++itr )
        os << Base::val2Char(*itr);
      return os; }

private:
    static unsigned const BITS_PER_BASE = 2;
    static storage_type const STORAGE_1 = storage_type(1);
    static storage_type const BASE_MASK =
            (STORAGE_1 << BITS_PER_BASE) - STORAGE_1;
    static unsigned const BITS_PER_STORAGE_UNIT =
            std::numeric_limits<storage_type>::digits;
    static unsigned const BASES_PER_STORAGE_UNIT =
            BITS_PER_STORAGE_UNIT / BITS_PER_BASE;
    static unsigned const STORAGE_UNITS_PER_KMER =
            (K+BASES_PER_STORAGE_UNIT-1)/BASES_PER_STORAGE_UNIT;
    static unsigned const UNUSED_TRAILING_BITS =
            BITS_PER_BASE*(BASES_PER_STORAGE_UNIT*STORAGE_UNITS_PER_KMER - K);
    static storage_type const FINAL_MASK =
            (STORAGE_1 << UNUSED_TRAILING_BITS) - STORAGE_1;
    static unsigned const SHIFT_LSBASE_TO_MSBASE =
            BITS_PER_BASE*(BASES_PER_STORAGE_UNIT-1u);

    void unused()
    {
        STATIC_ASSERT(K > 0u);
        STATIC_ASSERT(std::numeric_limits<storage_type>::is_integer);
        STATIC_ASSERT(!std::numeric_limits<storage_type>::is_signed);
        STATIC_ASSERT(BITS_PER_STORAGE_UNIT%BITS_PER_BASE == 0);
    }

    storage_type& data( size_type idx )
    { return mVal[idx/BASES_PER_STORAGE_UNIT]; }
    storage_type data( size_type idx ) const
    { return mVal[idx/BASES_PER_STORAGE_UNIT]; }
    static unsigned shift( size_type idx )
    { return BITS_PER_BASE *
                (BASES_PER_STORAGE_UNIT-1u-idx%BASES_PER_STORAGE_UNIT); }

    storage_type mVal[STORAGE_UNITS_PER_KMER];
};

template <unsigned K, class S>
struct Serializability< KMer<K,S> > : public TriviallySerializable {};

#endif /* KMER_H_ */
