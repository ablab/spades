///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file SmallVecDefs.h
 * \author tsharpe
 * \date Jul 8, 2009
 *
 * \brief Non-inline members for SmallVec class.
 */
#ifndef FEUDAL_SMALLVECDEFS_H_
#define FEUDAL_SMALLVECDEFS_H_

#include "feudal/SmallVec.h"
#include <cstddef>
#include <cstring>

template <class T, class A>
void SmallVec<T,A>::readFeudal( BinaryReader& reader,
                                unsigned long dataLen,
                                void* fixed )
{
    size_type size;
    // if the type T has a fixed external, binary-file size, then we can infer
    // our size from the dataLen arg.
    unsigned long elementSize = BinaryReader::externalSizeof(static_cast<T*>(0));
    if ( elementSize )
    {
        AssertEq( dataLen % elementSize, 0u );
        size = dataLen / elementSize;
    }
    // otherwise we'll have written our size into the fixed-data side stream
    else
    {
        Assert(fixed != 0);
        memcpy(&size,fixed,sizeof(size));
    }
    resize(size);
    reader.read(data(),dataEnd());
}

template <class T, class A>
void SmallVec<T,A>::init( T* last )
{
    T* first = dataEnd();
    mSize += last-first;
    A& alloc = allocator();
    while ( first != last )
        new (first++) T();
}

template <class T, class A>
void SmallVec<T,A>::init( T* last, T const& exemplar )
{
    T* first = dataEnd();
    mSize += last-first;
    A& alloc = allocator();
    while ( first != last )
        alloc.construct(first++,exemplar);
}

template <class T, class A>
void SmallVec<T,A>::destroy( T* first )
{
    T* last = dataEnd();
    mSize -= last-first;
    A& alloc = allocator();
    while ( last != first )
        alloc.destroy(--last);
}

template <class T, class A>
T* SmallVec<T,A>::swap( T* first, T* dest )
{
    using std::swap;
    T* last = dataEnd();
    while ( first != last )
        swap(*first++,*dest++);
    return dest;
}

template <class T, class A>
void SmallVec<T,A>::realloc( size_t nElements )
{
    ForceAssertLe(nElements,max_size());
    A& alloc = allocator();
    T* elements = nElements ? alloc.allocate(nElements,0) : 0;
    T* last = elements + size();
    T* src = dataEnd();
    while ( last != elements )
    {
        alloc.construct(--last,*--src);
        alloc.destroy(src);
    }
    deallocate();
    setData(elements);
    mCapacity = nElements;
}

// called to accomplish a swap when allocators are unequal
template <class T, class A>
void SmallVec<T,A>::exchange( SmallVec& that )
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
        swap(*src1++,*src2++);

    A& alloc = allocator();
    // copy and destroy the elements in this not matched by an element in that
    end = dataEnd();
    while ( src1 != end )
    {
        alloc.construct(src2++,*src1);
        alloc.destroy(src1++);
    }

    // copy and destroy the elements in that not matched by an element in this
    end = that.dataEnd();
    while ( src2 < end )
    {
        alloc.construct(src1++,*src2);
        alloc.destroy(src2++);
    }

    swap(mSize,that.mSize);
}

#endif /* FEUDAL_SMALLVECDEFS_H_ */
