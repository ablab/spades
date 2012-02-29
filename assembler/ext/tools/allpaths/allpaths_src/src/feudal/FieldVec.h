///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file FieldVec.h
 * \author ghall
 * \date Aug 11, 2009
 *
 * \brief A small (2^32 elements) vector that subdivides a byte to store a value.
 * This is an STL-compatible, vector-like container.
 */
#ifndef FEUDAL_FIELDVEC_H
#define FEUDAL_FIELDVEC_H

#include "Compare.h"
#include "CoreTools.h"
#include "feudal/BinaryStream.h"
#include "feudal/Generic.h"
#include "feudal/IsSizeT.h"
#include "feudal/Iterator.h"
#include "feudal/Oob.h"
#include "system/Assert.h"
#include "system/StaticAssert.h"
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <limits>
#include <cstring>

/// A vector whose elements are a subdivided byte.
/// Template args are N, the number of bits per value, and A, an allocator.
template <int N, class A>
class FieldVec
{
public:
    typedef unsigned char value_type;
    typedef value_type& reference;
    typedef value_type const& const_reference;
    typedef value_type* pointer;
    typedef value_type const* const_pointer;
    typedef unsigned int size_type;
    typedef std::ptrdiff_t difference_type;
    typedef std::iterator<std::random_access_iterator_tag,
                          value_type,
                          difference_type,
                          void,
                          value_type> ItrBase;
    typedef A allocator_type;

    class iterator
    : public ItrBase,
      public IteratorBase<iterator,size_type,difference_type>
    {
    public:
        iterator() : mpContainer(0) {}
        iterator( FieldVec* pContainer, size_type pos )
        : IteratorBase<iterator,size_type,difference_type>(pos),
          mpContainer(pContainer) {}

        // compiler-supplied copying and destructor are OK

        value_type operator*() const { return (*mpContainer)[this->pos()]; }

        value_type operator[]( difference_type idx ) const
        { return (*mpContainer)[this->pos()+idx]; }

        void set( value_type val )
        { mpContainer->set(this->pos(),val); }

    protected:
        FieldVec* mpContainer;
    };

    class const_iterator
    : public ItrBase,
      public IteratorBase<const_iterator,size_type,difference_type>
    {
    public:
        const_iterator() : mpContainer(0) {}
        const_iterator( FieldVec const* pContainer, size_type pos )
        : IteratorBase<const_iterator,size_type,difference_type>(pos),
          mpContainer(pContainer) {}

        // compiler-supplied copying and destructor are OK

        value_type operator*() const { return (*mpContainer)[this->pos()]; }
        value_type operator[]( difference_type idx ) const
        { return (*mpContainer)[this->pos()+idx]; }

    private:
        FieldVec const* mpContainer;
    };

    class lvalue
    {
    public:
        lvalue( FieldVec* pContainer, size_type pos )
        : mPos(pos), mpContainer(pContainer) {}

        value_type operator=( lvalue const& lv )
        { value_type val = lv; mpContainer->set(mPos,val); return val; }

        // compiler-supplied copying and destructor are OK

        operator value_type() const { return (*mpContainer)[mPos]; }

        value_type operator=( value_type val )
        { mpContainer->set(mPos,val); return val; }

        value_type operator+=( value_type val )
        { val += *this; mpContainer->set(mPos,val); return val; }

        // TODO: might consider other assignment operators as the need arises

    private:
        size_type mPos;
        FieldVec* mpContainer;
    };

    class lvalue_iterator
    : public ItrBase,
      public IteratorBase<lvalue_iterator,size_type,difference_type>
    {
    public:
        lvalue_iterator() : mpContainer(0) {}
        lvalue_iterator( FieldVec* pContainer, size_type pos )
        : IteratorBase<lvalue_iterator,size_type,difference_type>(pos),
          mpContainer(pContainer) {}

        // compiler-supplied copying and destructor are OK

        lvalue operator*() const { return lvalue(mpContainer,this->pos()); }

        lvalue operator[]( difference_type idx ) const
        { return lvalue(mpContainer,this->pos()+idx); }

    protected:
        FieldVec* mpContainer;
    };

    class reverse_iterator
    : public ItrBase,
      public IteratorBase<reverse_iterator,size_type,difference_type>
    {
    public:
        reverse_iterator() : mpContainer(0), mLast(~0) {}
        reverse_iterator( FieldVec* pContainer, size_type pos )
        : IteratorBase<reverse_iterator,size_type,difference_type>(pos),
          mpContainer(pContainer), mLast(pContainer->size()-1) {}

        // compiler-supplied copying and destructor are OK

        value_type operator*() const { return (*mpContainer)[mLast-this->pos()]; }

        value_type operator[]( difference_type idx ) const
        { return (*mpContainer)[mLast-(this->pos()+idx)]; }

        void set( value_type val )
        { mpContainer->set(mLast-this->pos(),val); }

    protected:
        FieldVec* mpContainer;
        size_type mLast;
    };

    class const_reverse_iterator
    : public ItrBase,
      public IteratorBase<const_reverse_iterator,size_type,difference_type>
    {
    public:
        const_reverse_iterator() : mpContainer(0), mLast(~0) {}
        const_reverse_iterator( FieldVec const* pContainer, size_type pos )
        : IteratorBase<const_reverse_iterator,size_type,difference_type>(pos),
          mpContainer(pContainer), mLast(pContainer->size()-1) {}

        // compiler-supplied copying and destructor are OK

        value_type operator*() const { return (*mpContainer)[mLast-this->pos()]; }

        value_type operator[]( difference_type idx ) const
        { return (*mpContainer)[mLast-(this->pos()+idx)]; }

    protected:
        FieldVec const* mpContainer;
        size_type mLast;
    };

    struct NopMapper
    { value_type operator()( value_type val ) const { return val; } };

    FieldVec()
    : mSize(0), mCapacity(0), mData(0)
    {
        allocInit(A());
    }

    explicit FieldVec( A const& alloc )
    : mSize(0), mCapacity(0), mData(0)
    {
        allocInit(alloc);
    }

    explicit FieldVec( size_type sz, value_type exemplar = 0,
                        size_type cap = 0, A const& alloc = A() )
    : mSize(0), mCapacity(0), mData(0)
    {
        allocInit(alloc);
        using std::max;
        reserve(max(sz, cap));
        init(sz, exemplar);
    }

    template <class Itr>
    FieldVec( Itr first, Itr const& last,
                size_type cap = 0, A const& alloc = A() )
    : mSize(0), mCapacity(0), mData(0)
    {
        allocInit(alloc);
        assign(first,last,cap);
    }

    template <class Itr, class Mapper>
    FieldVec( Itr first, Itr const& last, Mapper mapper,
              size_type cap = 0, A const& alloc = A() )
    : mSize(0), mCapacity(0), mData(0)
    {
        allocInit(alloc);
        using std::distance; size_t dist = distance(first,last);
        using std::max; reserve(max(dist,static_cast<size_t>(cap)));
        itercopy(first, last, mapper);
    }

    FieldVec( FieldVec const& that )
    : mSize(0), mCapacity(0), mData(0)
    {
        allocInit(A());
        reserve(that.size());
        bytecopy(that.data(), that.dataEnd());
        mSize = that.size();
    }

    ~FieldVec()
    {
        UtilizationReporter::gInstance.report(this,N,size(),capacity(),"FV");
        clear();
        deallocate();
        allocator().~A();
    }

    template <class A1>
    FieldVec& operator=( FieldVec<N,A1> const& that )
    {
        if ( this != &that )
        {
            clear();
            reserve(that.size());
            bytecopy(that.data(), that.dataEnd());
            mSize = that.size();
        }
        return *this;
    }

    FieldVec& operator=( FieldVec const& that )
    {
        if ( this != &that )
        {
            clear();
            reserve(that.size());
            bytecopy(that.data(), that.dataEnd());
            mSize = that.size();
        }
        return *this;
    }

    //
    // iterators
    //

    iterator begin() { return iterator(this,0); }
    iterator begin( size_type idx )
    { AssertLe(idx,size()); return iterator(this,idx); }
    iterator end() { return iterator(this,size()); }

    lvalue_iterator lbegin() { return lvalue_iterator(this,0); }
    lvalue_iterator lbegin( size_type idx )
    { AssertLe(idx,size()); return lvalue_iterator(this,idx); }
    lvalue_iterator lend() { return lvalue_iterator(this,size()); }

    const_iterator begin() const { return const_iterator(this,0); }
    const_iterator begin( size_type idx ) const
    { AssertLe(idx,size()); return const_iterator(this,idx); }
    const_iterator end() const { return const_iterator(this,size()); }

    const_iterator cbegin() const { return const_iterator(this,0); }
    const_iterator cbegin( size_type const idx ) const
    { AssertLe(idx,size()); return const_iterator(this,idx); }
    const_iterator cend() const { return const_iterator(this,size()); }

    reverse_iterator rbegin() { return reverse_iterator(this,0); }
    reverse_iterator rbegin( size_type idx )
    { AssertLe(idx,size()); return reverse_iterator(this,idx); }
    reverse_iterator rend() { return reverse_iterator(this, size()); }

    const_reverse_iterator rbegin() const
    { return const_reverse_iterator(this,0); }
    const_reverse_iterator rbegin( size_type idx ) const
    { AssertLe(idx,size()); return const_reverse_iterator(this,idx); }
    const_reverse_iterator rend() const
    { return const_reverse_iterator(this,size()); }

    const_reverse_iterator crbegin() const
    { return const_reverse_iterator(this,0); }
    const_reverse_iterator crbegin( size_type idx ) const
    { AssertLe(idx,size()); return const_reverse_iterator(this,idx); }
    const_reverse_iterator crend() const { return const_reverse_iterator(this, size()); }

    //
    // counts and sizes
    //

    FieldVec& resize(size_type sz, value_type exemplar = 0 )
    {
        reserve(sz);
        if (size() > sz)
        {
            mSize = sz;
        }
        else if (size() < sz)
        {
            init(sz, exemplar);
        }
        return *this;
    }

    size_type max_size() const
    {
        return std::numeric_limits<size_type>::max() - 1000;
    }

    bool empty() const
    {
        return !size();
    }

    bool full() const
    {
        return mSize == mCapacity;
    }

    size_type size() const
    {
        return mSize;
    }

    size_type allocSize() const
    {
        return physicalSize(mSize);
    }

    size_type capacity() const
    {
        return mCapacity;
    }

    FieldVec& reserve( size_type sz )
    {
        if ( sz > capacity() )
        {
            realloc(sz);
        }
        return *this;
    }

    //
    // element access
    //

    value_type operator[]( size_type idx ) const
    {
        return getValue(idx);
    }

    value_type at( size_type idx ) const
    {
        if ( idx >= size() )
        {
            OutOfBoundsReporter::oob("FieldVec", idx, size());
        }
        return getValue(idx);
    }

    value_type front() const
    {
        return getValue(0);
    }

    value_type back() const
    {
        return getValue(size() - 1);
    }

    //
    // modifiers
    //

    void set( size_type idx, value_type val )
    {
        if ( idx >= size() )
        {
            OutOfBoundsReporter::oob("FieldVec", idx, size());
        }
        Assert(!(val&~VAL_MASK));
        setValue(idx, val);
    }

    FieldVec& assign( size_type sz, value_type exemplar, size_type cap = 0 )
    {
        return assignConst(sz,exemplar,cap);
    }

    template <class Itr>
    FieldVec& assign( Itr first, Itr const& last, size_type cap = 0 )
    {
        return assignItr(first,last,cap,IsSizeT<Itr>());
    }

    template <class Itr, class Mapper>
    FieldVec& assign( Itr first, Itr const& last,
                        Mapper mapper, size_type cap = 0 )
    {
        clear();
        using std::distance; size_t dist = distance(first,last);
        using std::max; reserve(max(dist, static_cast<size_t>(cap)));
        itercopy(first, last, mapper);
        return *this;
    }

    FieldVec& append( FieldVec const& that )
    {
        return append(that.begin(),that.end());
    }

    FieldVec& append( size_type sz, value_type exemplar, size_type cap = 0 )
    {
        return appendConst(sz,exemplar,cap);
    }

    template <class Itr>
    FieldVec& append( Itr first, Itr const& last, size_type cap = 0 )
    {
        return appendItr(first,last,cap,IsSizeT<Itr>());
    }

    template <class Itr, class Mapper>
    FieldVec& append( Itr first, Itr const& last, Mapper mapper,
                        size_type cap = 0 )
    {
        using std::distance;
        using std::max;

        size_t newSize = size() + static_cast<size_t>(distance(first,last));
        reserve(max(newSize,static_cast<size_t>(cap)));
        itercopy(first, last, mapper);
        return *this;
    }

    FieldVec& push_back( value_type val )
    {
        return push_back(val, 1.6f, 64U);
    }

    FieldVec& push_back( value_type val,
                        float growthFact, size_type growthIncr )
    {
        Assert(!(val&~VAL_MASK));
        AssertGe(growthFact,0.f);
        if (full())
        {
            size_t sz = static_cast<size_t>(growthFact * size() + growthIncr);
            using std::max; using std::min;
            realloc(max(size()+1ul,min(sz,static_cast<size_t>(max_size()))));
        }
        setValue(mSize++, val);
        return *this;
    }

    iterator erase( iterator itr )
    {
        realign(itr.pos() + 1, itr.pos());
        --mSize;
        return itr;
    }

    iterator erase( iterator first, iterator last )
    {
        difference_type diff = last - first;
        AssertGe(diff,0);
        realign(last.pos(), first.pos());
        mSize -= diff;
        return first;
    }

    FieldVec& pop_back()
    {
        Assert(!empty());
        --mSize;
        return *this;
    }

    FieldVec& clear()
    {
        mSize = 0;
        return *this;
    }

    FieldVec& reverse()
    {
        for ( size_type lo = 0, hi = size(), end = hi/2; lo != end; ++lo )
        {
            value_type tmp = getValue(lo);
            set(lo,getValue(--hi));
            set(hi,tmp);
        }
        return *this;
    }

    FieldVec& shrink_to_fit()
    {
        if ( physicalSize(size()) != physicalSize(capacity()) )
        {
            realloc(size());
        }
        return *this;
    }

    FieldVec& swap( FieldVec& that )
    {
        if ( allocator() != that.allocator() )
        {
            // we have to do a slow swap
            exchange(that);
        } else
        {
            // we can do a quick swap -- just exchange data pointers
            unsigned long flippedBits = (mData ^ that.mData) & MASK;
            mData ^= flippedBits;
            that.mData ^= flippedBits;
            using std::swap;
            swap(mSize,that.mSize);
            swap(mCapacity,that.mCapacity);
        }
        return *this;
    }

    friend int compare( FieldVec const& v1, FieldVec const& v2 )
    {
        const_iterator itr1 = v1.begin();
        const_iterator itr2 = v2.begin();
        using std::min; const_iterator end = itr1 + min(v1.size(),v2.size());
        int result = 0;
        while ( !result && itr1 != end )
        { result = ::compare(*itr1,*itr2); ++itr1; ++itr2; }
        if ( !result ) result = ::compare(v1.size(),v2.size());
        return result;
    }

    A get_allocator() const
    {
        return allocator();
    }

    size_t writeFeudal( BinaryWriter& writer, void const** ppFixed ) const
    { *ppFixed = &mSize;
      return writer.write(data(),dataEnd()); }

    void readFeudal( BinaryReader& rdr, size_t varDataLen, void* pFixed )
    { size_type sz; memcpy(&sz,pFixed,sizeof(sz));
      ForceAssertEq(physicalSize(sz),varDataLen);
      resize(sz);
      rdr.read(data(),dataEnd()); }

    size_t writeBinary( BinaryWriter& writer ) const
    { size_t len = writer.write(mSize);
      return len+writer.write(data(),dataEnd()); }

    void readBinary( BinaryReader& reader )
    { size_type sz; reader.read(&sz); resize(sz);
      reader.read(data(),dataEnd()); }

    static size_t externalSizeof() { return 0; }

    static unsigned int fixedDataLen()
    { return sizeof(size_type); }

    static size_type physicalSize( size_type sz )
    { return (sz + VALS_PER_BYTE - 1) / VALS_PER_BYTE; }

protected:

    static size_type logicalSize( size_type sz )
    { return sz * VALS_PER_BYTE; }

    pointer data()
    { return reinterpret_cast<pointer>(mData&MASK); }

    const_pointer data() const
    { return reinterpret_cast<const_pointer>(mData&MASK); }

    pointer dataEnd()
    { return data() + physicalSize(size()); }

    const_pointer dataEnd() const
    { return data() + physicalSize(size()); }

    void setSize( size_type sz )
    { Assert(sz <= mCapacity); mSize = sz; }

private:
    void setData( void const* data )
    {
        mData ^= (mData ^ reinterpret_cast<unsigned long>(data)) & MASK;
    }

    void allocInit( A const& alloc )
    {
        STATIC_ASSERT(sizeof(A) <= 2u);
        new (&allocator()) A(alloc);
    }

    A& allocator()
    {
        return reinterpret_cast<A*>(&mData + 1)[-1];
    }

    A const& allocator() const
    {
        return reinterpret_cast<A const*>(&mData + 1)[-1];
    }

    void init( size_type last, value_type exemplar );

    template <class Itr, class Mapper>
    void itercopy( Itr first, Itr const& last, Mapper mapper )
    {
        size_type offset = size();
        mSize += std::distance(first,last);
        while ( first != last )
        {
            value_type val = mapper(*first);
            Assert(!(val&~VAL_MASK));
            setValue(offset++, val);
            ++first;
        }
    }

    void bytecopy( const_pointer first, const_pointer last )
    {
        memcpy(dataEnd(), first, last - first);
    }

    FieldVec& assignConst( size_type sz, value_type exemplar, size_type cap )
    {
        using std::max;
        reserve(max(sz, cap));
        clear();
        init(sz, exemplar);
        return *this;
    }

    template <class Itr>
    FieldVec& assignItr( Itr first, Itr const& last, size_type cap, NoSizeT )
    {
        clear();
        using std::distance;
        size_t dist = distance(first,last);
        using std::max;
        reserve(max(dist, static_cast<size_t>(cap)));
        itercopy(first, last, NopMapper());
        return *this;
    }

    template <class Itr>
    FieldVec& assignItr( Itr first, Itr const& last, size_type cap, YesSizeT )
    {
        return assignConst(first,last,cap);
    }

    FieldVec& appendConst( size_type sz, value_type exemplar, size_type cap )
    {
        using std::max;
        reserve(max(size()+static_cast<size_t>(sz), static_cast<size_t>(cap)));
        init(size() + sz, exemplar);
        return *this;
    }

    template <class Itr>
    FieldVec& appendItr( Itr first, Itr const& last, size_type cap, NoSizeT )
    {
        using std::distance;
        using std::max;

        size_t newSize = size() + static_cast<size_t>(distance(first, last));
        reserve(max(newSize, static_cast<size_t>(cap)));
        itercopy(first, last, NopMapper());
        return *this;
    }

    template <class Itr>
    FieldVec& appendItr( Itr first, Itr const& last, size_type cap, YesSizeT )
    {
        return appendConst(first,last,cap);
    }

    void exchange( FieldVec& that );
    void realloc( size_t sz );

    void setValue( size_type idx, value_type val )
    {
        AssertLt( idx, size( ) );
        value_type& datum = data()[idx/VALS_PER_BYTE];
        size_type shift = idx % VALS_PER_BYTE * N;
        datum ^= (datum ^ (val << shift)) & (VAL_MASK << shift);
    }

    value_type getValue( size_type idx ) const
    {
        AssertLt( idx, size( ) );
        return (data()[idx/VALS_PER_BYTE] >> idx%VALS_PER_BYTE*N) & VAL_MASK;
    }

    void realign( size_type src, size_type dest );

    void deallocate()
    {
        if ( data() )
        {
            allocator().deallocate(data(), physicalSize(capacity()));
        }
    }

    void checkN();

    size_type mSize;
    size_type mCapacity;
    unsigned long mData;
    static unsigned long const MASK = 0xffffffffffff; // lowest 48 bits
    static size_type const VALS_PER_BYTE = 8/N;
    static value_type const VAL_MASK = (1U << N) - 1U;
};

template <int N, class A>
struct Serializability<FieldVec<N,A> > : public SelfSerializable {};

template <int N, class A>
bool operator==( FieldVec<N,A> const& v1, FieldVec<N,A> const& v2 )
{
    using std::equal;
    return ((v1.size() == v2.size()) &&
                (equal(v1.begin(), v1.end(), v2.begin())));
}


template <int N, class A1, class A2>
bool operator==( FieldVec<N,A1> const& v1, FieldVec<N,A2> const& v2 )
{
    using std::equal;
    return ((v1.size() == v2.size()) &&
                (equal(v1.begin(), v1.end(), v2.begin())));
}

template <int N, class A1, class A2>
bool operator!=( FieldVec<N,A1> const& v1, FieldVec<N,A2> const& v2 )
{
    using std::equal;
    return ((v1.size() != v2.size()) ||
                (!equal(v1.begin(), v1.end(), v2.begin())));
}

template <int N, class A1, class A2>
bool operator<( FieldVec<N,A1> const& v1, FieldVec<N,A2> const& v2 )
{
    using std::lexicographical_compare;
    return lexicographical_compare(v1.begin(), v1.end(), v2.begin(), v2.end());
}

template <int N, class A1, class A2>
bool operator>( FieldVec<N,A1> const& v1, FieldVec<N,A2> const& v2 )
{
    return (v2 < v1);
}


template <int N, class A1, class A2>
bool operator<=( FieldVec<N,A1> const& v1, FieldVec<N,A2> const& v2 )
{
    return !(v2 < v1);
}

template <int N, class A1, class A2>
bool operator>=( FieldVec<N,A1> const& v1, FieldVec<N,A2> const& v2 )
{
    return !(v1 < v2);
}

template <int N, class A>
void swap( FieldVec<N,A>& v1, FieldVec<N,A>& v2 )
{
    v1.swap(v2);
}

#endif  // FEUDAL_FIELDVEC_H
