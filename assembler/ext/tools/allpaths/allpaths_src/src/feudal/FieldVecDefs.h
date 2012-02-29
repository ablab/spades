///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef _FIELDVECDEFS_H_
#define _FIELDVECDEFS_H_

#include "feudal/FieldVec.h"
#include <cstddef>

template <int N, class A>
void FieldVec<N,A>::init( size_type last,
                          value_type exemplar )
{
    Assert(!(exemplar&~VAL_MASK));

    size_type first = size();
    mSize = last;

    // march up to the first byte boundary
    while ( first%VALS_PER_BYTE && first != last )
        setValue(first++, exemplar);

    // if there's more to do after filling a ragged initial byte
    if ( first != last )
    {
        // replicate the exemplar bit pattern throughout the byte
        if ( N == 1 )
            exemplar |= exemplar << 1;
        if ( N <= 2 )
            exemplar |= exemplar << 2;
        if ( N <= 4 )
            exemplar |= exemplar << 4;

        // fill the remainder of the bytes with that replicated pattern
        memset(data()+first/VALS_PER_BYTE,exemplar,physicalSize(last-first));
    }
}

template <int N, class A>
void FieldVec<N,A>::realloc( size_t nElements )
{
    ForceAssertLe(nElements,max_size());
    A& alloc = allocator();
    pointer elements = 0;
    size_type newSize = physicalSize(nElements);
    if ( newSize )
    {
        elements = alloc.allocate(newSize, 0);
        size_type curSize = physicalSize(size());
        if( curSize )
        {
            memcpy(elements, data(), curSize);
        }
    }
    deallocate();
    setData(elements);
    mCapacity = logicalSize(newSize);
}

template <int N, class A>
void FieldVec<N,A>::realign( size_type src,
                             size_type dest )
{
    size_type end = size();
    while ( src != end )
    {
        setValue(dest++, getValue(src++));
    }
}

// called to accomplish a swap when allocators are unequal
template <int N, class A>
void FieldVec<N,A>::exchange( FieldVec& that )
{
    using std::swap;

    // make sure there's adequate space in the shorter vector
    reserve(that.size());
    that.reserve(size());

    pointer src1 = data();
    pointer src2 = that.data();
    pointer end = src1 + physicalSize(std::min(size(),that.size()));

    // swap the places where we both have bytes
    while ( src1 != end )
    {
        swap(*src1++,*src2++);
    }

    // copy the bytes in this not matched by bytes in that
    size_type len;
    if ( (len = dataEnd() - src1) )
    {
        memcpy(src2,src1,len);
    }
    // copy the bytes in that not matched by bytes in this
    else if ( (len = that.dataEnd() - src2) )
    {
        memcpy(src1,src2,len);
    }

    // swap the sizes
    swap(mSize, that.mSize);
}

template <int N, class A>
void FieldVec<N,A>::checkN()
{
    STATIC_ASSERT( N==1 || N==2 || N==4 );
}

#endif /* _FIELDVECDEFS_H_ */
