///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file OuterVec.h
 * \author tsharpe
 * \date Jun 8, 2009
 *
 * \brief A vector that manages a sub-allocator on behalf of its elements.
 * This is an STL-compatible, vector-like container that manages memory for its
 * vector-like elements. The purpose is to avoid the overhead of a zillion tiny
 * memory allocations by the elements.
 * It's very similar to std::vector, with one important exception:  reallocation
 * happens by swapping elements from their old location to a new one (rather than
 * by copying, as with std::vector).  This means that the memory pool needn't be
 * expanded as a result of reallocation.
 */
#ifndef FEUDAL_OUTERVEC_H_
#define FEUDAL_OUTERVEC_H_

#include "feudal/Mempool.h"
#include "feudal/Iterator.h"
#include "feudal/Oob.h"
#include "feudal/BinaryStream.h"
#include "feudal/IsSizeT.h"
#include "feudal/Generic.h"
#include "system/Assert.h"
#include "Compare.h"
#include <algorithm>
#include <cstddef>
#include <limits>
#include <memory>

#if 0
// Here is the interface for T, the OuterVec's value_type.
class T
{
public:
    typedef anything value_type;
    typedef anyint size_type;
    T();
    T( T const& );
    T( A const& allocator );
    T& operator=( T const& );
    size_type size() const;
    void swap( T& that ); // must be a swap-in-place
    size_type allocSize() const;
};

// class S, the sub-allocator, must be default-constructable, it must look like
// a standard allocator, but must also have two additional methods:
// preAllocate(size_t,size_t) and getMaxEnchunkableSize().
// it must also be copy-constructable if you use the OuterVec constructor that
//   takes an S

#endif

template <class T,
          class S = MempoolOwner<typename T::value_type>,
          class A = std::allocator<T> >
class OuterVec
{
public:
    typedef T value_type;
    typedef T& reference;
    typedef T const& const_reference;
    typedef T* pointer;
    typedef T const* const_pointer;
    typedef size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef FwdIterator<T> iterator;
    typedef FwdConstIterator<T> const_iterator;
    typedef RevIterator<T> reverse_iterator;
    typedef RevConstIterator<T> const_reverse_iterator;

    OuterVec() : mpElements(0), mpCurEnd(0), mpAllocEnd(0) {}

    explicit OuterVec( S const& subAlloc )
    : mSubAllocator(subAlloc), mpElements(0), mpCurEnd(0), mpAllocEnd(0) {}

    explicit OuterVec( size_type sz )
    : mpElements(0), mpCurEnd(0), mpAllocEnd(0)
    { reserve(sz);
      mpCurEnd = init(data()+sz); }

    explicit OuterVec( size_type sz, T const& exemplar,
                        size_type cap = 0, A const& alloc = A() )
    : mAllocator(alloc), mpElements(0), mpCurEnd(0), mpAllocEnd(0)
    { using std::max; reserve(max(sz,cap));
      mpCurEnd = init(data()+sz,exemplar); }

    template <class Itr>
    OuterVec( Itr first, Itr const& last,
               size_type cap = 0, A const& alloc = A() )
    : mAllocator(alloc), mpElements(0), mpCurEnd(0), mpAllocEnd(0)
    { assign(first,last,cap); }

    OuterVec( OuterVec const& that )
    : mpElements(0), mpCurEnd(0), mpAllocEnd(0)
    { reserve(that.size());
      mpCurEnd = copy(that.data(),that.dataEnd()); }

    template <class S1, class A1>
    OuterVec( OuterVec<T,S1,A1> const& that )
    : mpElements(0), mpCurEnd(0), mpAllocEnd(0)
    { reserve(that.size());
      mpCurEnd = copy(that.begin(),that.end()); }

    ~OuterVec()
    { UtilizationReporter::gInstance.report(this,8*sizeof(T),size(),capacity(),"OV");
      clear();
      deallocate(); }

    template <class S1,class A1>
    OuterVec& operator=( OuterVec<T,S1,A1> const& that )
    { clear(); reserve(that.size());
      mpCurEnd = copy(that.begin(),that.end());
      return *this; }

    OuterVec& operator=( OuterVec const& that )
    { if ( this != &that )
      { clear(); reserve(that.size());
        mpCurEnd = copy(that.data(),that.dataEnd()); }
      return *this; }

    /*
     * iterators
     */

    iterator begin() { return iterator(data()); }
    iterator begin( size_type idx )
    { AssertLe(idx,size()); return iterator(data()+idx); }

    const_iterator begin() const { return const_iterator(data()); }
    const_iterator begin( size_type idx ) const
    { AssertLe(idx,size()); return const_iterator(data()+idx); }

    const_iterator cbegin() const { return const_iterator(data()); }
    const_iterator cbegin( size_type idx ) const
    { AssertLe(idx,size()); return const_iterator(data()+idx); }

    iterator end() { return iterator(dataEnd()); }
    const_iterator end() const { return const_iterator(dataEnd()); }
    const_iterator cend() const { return const_iterator(dataEnd()); }

    reverse_iterator rbegin() { return reverse_iterator(0,dataEnd()); }
    reverse_iterator rbegin( size_type idx )
    { return reverse_iterator(idx,dataEnd()); }

    const_reverse_iterator rbegin() const
    { return const_reverse_iterator(0,dataEnd()); }
    const_reverse_iterator rbegin( size_type idx ) const
    { AssertLe(idx,size()); return const_reverse_iterator(idx,dataEnd()); }

    const_reverse_iterator crbegin() const
    { return const_reverse_iterator(0,dataEnd()); }
    const_reverse_iterator crbegin( size_type idx ) const
    { return const_reverse_iterator(idx,dataEnd()); }

    reverse_iterator rend() { return reverse_iterator(size(),dataEnd()); }
    const_reverse_iterator rend() const
    { return const_reverse_iterator(size(),dataEnd()); }
    const_reverse_iterator crend() const
    { return const_reverse_iterator(size(),dataEnd()); }

    /*
     * counts and sizes
     */

    size_type size() const { return mpCurEnd - mpElements; }

    size_type allocSize() const { return size(); }

    /// sum of size() over each element
    size_t sumSizes() const
    { size_t sum = 0;
      const_iterator stop(end());
      for ( const_iterator itr(begin()); itr != stop; ++itr )
            sum += itr->size();
      return sum; }


    OuterVec& resize( size_type sz )
    { reserve(sz);
      T* posn = data() + sz;
      T* end = dataEnd();
      if ( end > posn ) mpCurEnd = destroy(posn);
      if ( posn > end ) mpCurEnd = init(posn);
      return *this; }

    OuterVec& resize( size_type sz, T const& exemplar )
    { reserve(sz);
      T* posn = data() + sz;
      T* end = dataEnd();
      if ( end > posn ) mpCurEnd = destroy(posn);
      if ( posn > end ) mpCurEnd = init(posn,exemplar);
      return *this; }

    size_type max_size() const
    { return 1UL << 48; }

    size_type capacity() const { return mpAllocEnd - mpElements; }
    bool empty() const { return mpCurEnd == mpElements; }
    bool full() const { return mpCurEnd == mpAllocEnd; }

    OuterVec& reserve( size_type sz )
    { if ( sz > capacity() ) realloc(sz); return *this; }

    /*
     * element access
     */

    reference operator[]( size_type idx )
    { AssertLt(idx, size()); return data()[idx]; }

    const_reference operator[]( size_type idx ) const
    { AssertLt(idx,size()); return data()[idx]; }

    reference at( size_type idx )
    { if ( idx >= size() ) OutOfBoundsReporter::oob("OuterVec",idx,size());
      return data()[idx]; }

    const_reference at( size_type idx ) const
    { if ( idx >= size() ) OutOfBoundsReporter::oob("OuterVec",idx,size());
      return data()[idx]; }

    reference front() { Assert(!empty()); return data()[0]; }
    const_reference front() const { Assert(!empty()); return data()[0]; }
    reference back() { Assert(!empty()); return dataEnd()[-1]; }
    const_reference back() const { Assert(!empty()); return dataEnd()[-1]; }

    /*
     * modifiers
     */

    OuterVec& assign( size_type sz )
    { return clear().resize(sz); }

    OuterVec& assign( size_type sz, T const& exemplar, size_type cap = 0 )
    { return assignConst(sz,exemplar,cap); }

    template <class Itr>
    OuterVec& assign( Itr first, Itr const& last, size_type cap = 0 )
    { return assignItr(first,last,cap,IsSizeT<Itr>()); }

    OuterVec& append( size_type sz )
    { resize(size()+sz);
      return *this; }

    OuterVec& append( size_type sz,
                      T const& exemplar,
                      size_type cap = 0 )
    { using std::max; size_type reserveSize = max(size()+sz,cap);
      if ( reserveSize > capacity() && isSelf(exemplar) )
      { append(sz,T(exemplar),cap); }
      else
      { reserve(reserveSize);
        mpCurEnd = init(dataEnd()+sz,exemplar); }
      return *this; }

    template <class Itr>
    OuterVec& append( Itr first, Itr const& last, size_type cap = 0 )
    { Assert(first==last||!isSelf(*first));
      using std::max; using std::distance;
      reserve(max(size()+static_cast<size_type>(distance(first,last)),cap));
      mpCurEnd = copy(first,last);
      return *this; }

    OuterVec& push_back( T const& val )
    { return push_back(val,1.6f,64UL); }

    OuterVec& push_back( T const& val, float growthFact, size_type growthIncr )
    { if ( !full() )
      { mpCurEnd = init(dataEnd()+1,val); }
      else if ( isSelf(val) )
      { push_back(T(val),growthFact,growthIncr); }
      else
      { size_type sz = static_cast<size_type>(growthFact*size()+growthIncr);
        using std::max; realloc(max(size()+1,sz));
        mpCurEnd = init(dataEnd()+1,val); }
      return *this; }

    iterator erase( iterator pos )
    { T* dest = &*pos;
      Assert(data() <= dest && dest <= dataEnd());
      mpCurEnd = destroy(swap(dest+1,dest));
      return pos; }

    iterator erase( iterator first, iterator const& last )
    { T* dest = &*first;
      T* start = &*last;
      Assert(data() <= dest && dest <= start && start <= dataEnd());
      mpCurEnd = destroy(swap(start,dest));
      return first; }

    /// Erases each element for which *Itr is true.
    template <class Itr>
    OuterVec& eraseIf( Itr begin, Itr const& end )
    { using std::distance; using std::swap;
      AssertEq(static_cast<size_type>(distance(begin,end)),size());
      T* dest = data();
      for ( T* src = dest; begin != end; ++begin, ++src )
      { if ( !*begin ) { if ( src != dest ) swap(*src,*dest); ++dest; } }
      mpCurEnd = destroy(dest);
      return *this; }

    /// Erases each element for which *Itr is false.
    template <class Itr>
    OuterVec& eraseUnless( Itr begin, Itr const& end )
    { using std::distance; using std::swap;
      AssertEq(static_cast<size_type>(distance(begin,end)),size());
      T* dest = data();
      for ( T* src = dest; begin != end; ++begin, ++src )
      { if ( *begin ) { if ( src != dest ) swap(*src,*dest); ++dest; } }
      mpCurEnd = destroy(dest);
      return *this; }

    /// Erases each element with an index called out by *Itr.
    /// Indices to erase must be in strictly increasing order.
    template <class Itr>
    OuterVec& eraseEntries( Itr begin, Itr const& end )
    { using std::swap;
      T* elements = data();
      size_type nnn = size();
      size_type dest = 0;
      for ( size_type iii = 0; iii < nnn; ++iii )
      { if ( begin != end && static_cast<size_type>(*begin) == iii )
        { ++begin; Assert(begin==end||static_cast<size_type>(*begin)>iii); }
        else
        {  if ( iii != dest ) swap(elements[iii],elements[dest]); ++dest; } }
      mpCurEnd = destroy(elements+dest);
      return *this; }

    OuterVec& pop_back()
    { Assert(!empty());
      mpCurEnd = destroy(dataEnd()-1);
      return *this; }

    OuterVec& clear() { mpCurEnd = destroy(data()); return *this; }
    OuterVec& shrink_to_fit() { if ( !full() ) realloc(size()); return *this; }

    OuterVec& swap( OuterVec& that )
    { // if we have to do a slow swap because allocators are unequal
      if ( mAllocator != that.mAllocator )
        exchange(that);
      else // we can do a quick swap -- just exchange data pointers
      { using std::swap;
        swap(mSubAllocator,that.mSubAllocator);
        swap(mpElements,that.mpElements);
        swap(mpCurEnd,that.mpCurEnd);
        swap(mpAllocEnd,that.mpAllocEnd); }
      return *this; }

    friend int compare( OuterVec const& v1, OuterVec const& v2 )
    { value_type* itr1 = v1.data();
      value_type* itr2 = v2.data();
      using std::min; value_type* end = itr1 + min(v1.size(),v2.size());
      int result = 0;
      while ( !result && itr1 != end ) result = ::compare(*itr1++,*itr2++);
      if ( !result ) result = ::compare(v1.size(),v2.size());
      return result; }

    A get_allocator() const { return mAllocator; }

    size_t writeBinary( BinaryWriter& writer ) const
    { size_t len = writer.write(size());
      return len+writer.write(data(),dataEnd()); }

    void readBinary( BinaryReader& reader )
    { size_type sz; reader.read(&sz); resize(sz,T(mSubAllocator));
      reader.read(data(),dataEnd()); }

    static size_t externalSizeof() { return 0; }

    size_t writeFeudal( BinaryWriter& writer, void const** ppFixed ) const
    { *const_cast<size_t*>(static_cast<size_t const*>(*ppFixed)) = size();
      return writer.write(data(),dataEnd()); }

    void readFeudal( BinaryReader& rdr, unsigned long dataLen, void* pFixed );

    static unsigned int fixedDataLen()
    { return BinaryReader::externalSizeof(static_cast<T*>(0)) ?
                0U : sizeof(size_type); }

protected:
    S& getSubAllocator()
    { return mSubAllocator; }
    S const& getSubAllocator() const
    { return mSubAllocator; }

private:
    T* data() { return mpElements; }
    T const* data() const { return mpElements; }
    T* dataEnd() { return mpCurEnd; }
    T const* dataEnd() const { return mpCurEnd; }

    bool isSelf( T const& val ) { return data() <= &val && &val < dataEnd(); }

    T* init( T* last );
    T* init( T* last, T const& exemplar );
    T* destroy( T* first );
    T* swap( T* first, T* dest );

    template <class Itr>
    T* copy( Itr first, Itr const& last )
    { T* dest = dataEnd();
      size_t nVs = 0;
      size_t maxVs = mSubAllocator.getMaxEnchunkableSize();
      for ( Itr itr(first); itr != last; ++itr )
      { size_t nnn = itr->allocSize(); if ( nnn <= maxVs ) nVs += nnn; }
      mSubAllocator.preAllocate(1UL,nVs);
      while ( first != last )
      { *(new (dest++) T(mSubAllocator)) = *first;
        ++first; }
      return dest; }

    OuterVec& assignConst( size_type sz, T const& exemplar, size_type cap )
    { using std::max; size_type reserveSize = max(sz,cap);
      if ( reserveSize > capacity() && isSelf(exemplar) )
      { assign(sz,T(exemplar),cap); }
      else
      { clear();
        reserve(reserveSize);
        mpCurEnd = init(data()+sz,exemplar); }
      return *this; }

    template <class Itr>
    OuterVec& assignItr( Itr first, Itr const& last, size_type cap, NoSizeT )
    { Assert(first==last||!isSelf(*first));
      clear();
      using std::max; using std::distance;
      reserve(max(static_cast<size_type>(distance(first,last)),cap));
      mpCurEnd = copy(first,last);
      return *this; }

    template <class Itr>
    OuterVec& assignItr( Itr first, Itr const& last, size_type cap, YesSizeT )
    { return assignConst(first,last,cap); }

    OuterVec& appendConst( size_type sz, T const& exemplar, size_type cap )
    { using std::max; size_type reserveSize = max(size()+sz,cap);
      if ( reserveSize > capacity() && isSelf(exemplar) )
      { append(sz,T(exemplar),cap); }
      else
      { reserve(reserveSize);
        mpCurEnd = init(dataEnd()+sz,exemplar); }
      return *this; }

    template <class Itr>
    OuterVec& appendItr( Itr first, Itr const& last, size_type cap, NoSizeT )
    { Assert(first==last||!isSelf(*first));
      using std::max; using std::distance;
      reserve(max(size()+static_cast<size_type>(distance(first,last)),cap));
      mpCurEnd = copy(first,last);
      return *this; }

    template <class Itr>
    OuterVec& appendItr( Itr first, Itr const& last, size_type cap, YesSizeT )
    { return appendConst(first,last,cap); }

    void realloc( size_type sz );
    void exchange( OuterVec& that );

    void deallocate()
    { if ( data() ) mAllocator.deallocate(data(),capacity()); }

    S mSubAllocator;
    A mAllocator;
    T* mpElements;
    T* mpCurEnd;
    T* mpAllocEnd;
};

template <class T, class S, class A>
struct Serializability< OuterVec<T,S,A> > : public SelfSerializable {};

template <class T, class S, class A1, class A2>
bool operator==( OuterVec<T,S,A1> const& v1, OuterVec<T,S,A2> const& v2 )
{ using std::equal;
  return v1.size()==v2.size() && equal(v1.begin(),v1.end(),v2.begin()); }

template <class T, class S, class A1, class A2>
bool operator!=( OuterVec<T,S,A1> const& v1, OuterVec<T,S,A2> const& v2 )
{ using std::equal;
  return v1.size()!=v2.size() || !equal(v1.begin(),v1.end(),v2.begin()); }

template <class T, class S, class A1, class A2>
bool operator<( OuterVec<T,S,A1> const& v1, OuterVec<T,S,A2> const& v2 )
{ using std::lexicographical_compare;
  return lexicographical_compare(v1.begin(),v1.end(),v2.begin(),v2.end()); }

template <class T, class S, class A1, class A2>
bool operator>( OuterVec<T,S,A1> const& v1, OuterVec<T,S,A2> const& v2 )
{ return (v2 < v1); }

template <class T, class S, class A1, class A2>
bool operator<=( OuterVec<T,S,A1> const& v1, OuterVec<T,S,A2> const& v2 )
{ return !(v2 < v1); }

template <class T, class S, class A1, class A2>
bool operator>=( OuterVec<T,S,A1> const& v1, OuterVec<T,S,A2> const& v2 )
{ return !(v1 < v2); }

template <class T, class S, class A>
void swap( OuterVec<T,S,A>& v1, OuterVec<T,S,A>& v2 )
{ v1.swap(v2); }

#endif /* FEUDAL_OUTERVEC_H_ */
