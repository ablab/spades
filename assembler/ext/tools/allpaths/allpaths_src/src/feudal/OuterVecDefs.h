///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file OutervecDefs.h
 * \author tsharpe
 * \date Jul 2, 2009
 *
 * \brief Definitions for non-inline OuterVec methods.
 */
#ifndef FEUDAL_OUTERVECDEFS_H_
#define FEUDAL_OUTERVECDEFS_H_

#include "feudal/OuterVec.h"
#include <cstddef>

template <class T, class S, class A>
void OuterVec<T,S,A>::readFeudal( BinaryReader& rdr, unsigned long dataLen,
                                    void* pFixed )
{
    size_type size;
    // if the type T has a fixed external, binary-file size, then we can infer
    // our size from the dataLen arg.
    size_t elementSize = BinaryReader::externalSizeof(static_cast<T*>(0));
    if ( elementSize )
    {
        AssertEq( dataLen % elementSize, 0u );
        size = dataLen / elementSize;
    }
    // otherwise we'll have written our size into the fixed-data side stream
    else
    {
        Assert(pFixed != 0);
        memcpy(&size,pFixed,sizeof(size));
    }
    resize(size,T(mSubAllocator));
    rdr.read(data(),dataEnd());
}

template <class T, class S, class A>
T* OuterVec<T,S,A>::init( T* last )
{
    T* first = dataEnd();
    while ( first != last )
        new (first++) T(mSubAllocator);
    return last;
}

template <class T, class S, class A>
T* OuterVec<T,S,A>::init( T* last, T const& exemplar )
{
    T* first = dataEnd();
    mSubAllocator.preAllocate(exemplar.allocSize(),last-first);
    while ( first != last )
        *(new (first++) T(mSubAllocator)) = exemplar;
    return last;
}

template <class T, class S, class A>
T* OuterVec<T,S,A>::destroy( T* first )
{
    T* last = dataEnd();
    while ( last != first )
        mAllocator.destroy(--last);
    return first;
}

template <class T, class S, class A>
T* OuterVec<T,S,A>::swap( T* first, T* dest )
{
    using std::swap;
    T* last = dataEnd();
    while ( first != last )
        swap(*first++,*dest++);
    return dest;
}

template <class T, class S, class A>
void OuterVec<T,S,A>::realloc( OuterVec<T,S,A>::size_type nElements )
{
    ForceAssertLe(nElements,max_size());
    T* elements = nElements ? mAllocator.allocate(nElements,0) : 0;
    T* last = elements + size();
    T* src = dataEnd();

    mpCurEnd = last;

    while ( last != elements )
    {
        (new (--last) T(mSubAllocator))->swap(*--src);
        mAllocator.destroy(src);
    }

    deallocate();
    mpElements = elements;
    mpAllocEnd = elements + nElements;
}

// called to accomplish a swap when allocators are unequal
template <class T, class S, class A>
void OuterVec<T,S,A>::exchange( OuterVec<T,S,A>& that )
{
    using std::swap;

    // make sure there's adequate space in the shorter vector
    reserve(that.size());
    that.reserve(size());

    T* src1 = data();
    T* src2 = that.data();
    T* end = src1 + std::min(size(),that.size());

    // swap the places where we both have elements
    while ( src1 != end )
        src1++->swap(*src2++);

    // copy and destroy the elements in this not matched by an element in that
    end = dataEnd();
    while ( src1 != end )
    {
        (new (src2++) T(mSubAllocator))->swap(*src1);
        mAllocator.destroy(src1++);
    }

    // copy and destroy the elements in that not matched by an element in this
    end = that.dataEnd();
    while ( src2 < end )
    {
        (new (src1++) T(mSubAllocator))->swap(*src2);
        mAllocator.destroy(src2++);
    }

    size_type siz = that.size();
    that.mpCurEnd = that.data() + size();
    mpCurEnd = data() + siz;
}

#define INSTANTIATE_OUTERVEC_INTERNALS(T,S) \
template T* OuterVec<T,S,std::allocator<T> >::init( T* last ); \
template T* OuterVec<T,S,std::allocator<T> >::init( T* last, T const& exemplar ); \
template T* OuterVec<T,S,std::allocator<T> >::destroy( T* first ); \
template T* OuterVec<T,S,std::allocator<T> >::swap( T* first, T* dest ); \
template void OuterVec<T,S,std::allocator<T> >::realloc( OuterVec<T,S,std::allocator<T> >::size_type nElements ); \
template void OuterVec<T,S,std::allocator<T> >::exchange( OuterVec<T,S,std::allocator<T> >& that )

#endif /* FEUDAL_OUTERVECDEFS_H_ */
