///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file HugeBVec.h
 * \author tsharpe
 * \date Sep 16, 2010
 *
 * \brief
 */
#ifndef FEUDAL_HUGEBVEC_H_
#define FEUDAL_HUGEBVEC_H_

#include "feudal/Iterator.h"
#include "system/Assert.h"
#include <cstddef>
#include <cstring>
#include <fstream>

/// A bvec that can be as big as all memory.
/// It's effectively read-only.
/// It's created by writing into a file using a builder sub-class, and then the
/// file is memory-mapped for use.
class HugeBVec
{
public:
    typedef unsigned char value_type;
    typedef size_t size_type;
    typedef std::ptrdiff_t difference_type;

    class Builder
    {
    public:
        Builder( char const* fileName )
        : mOS(fileName,
                std::ios_base::out|std::ios_base::binary|std::ios_base::trunc),
          mSize(0), mByte(0)
        { writeSize(); }

        ~Builder() { if ( mOS.is_open() ) close(); }

        void close()
        { if ( mSize & 3 ) { mByte <<= 2*(4-(mSize&3)); mOS << mByte; }
          mOS.seekp(0); writeSize(); mOS.close(); }

        size_type size() const { return mSize; }

        void push_back( value_type baseCode )
        { AssertLt(baseCode,4); mByte <<= 2; mByte |= baseCode&3;
          if ( !(++mSize & 3) ) mOS << mByte; }

        template <class Itr>
        void append( Itr itr, Itr const& end )
        { while ( itr != end ) { push_back(*itr); ++itr; } }

    private:
        Builder( Builder const& ); // unimplemented -- no copying
        Builder& operator=( Builder const& ); // unimplemented -- no copying

        void writeSize()
        { mOS.write(reinterpret_cast<char const*>(&mSize),sizeof(mSize)); }

        std::ofstream mOS;
        size_type mSize;
        value_type mByte;
    };

    typedef std::iterator<std::random_access_iterator_tag,
                          value_type,
                          difference_type,
                          void,
                          value_type> ItrTagBase;
    class const_iterator
    : public ItrTagBase,
      public IteratorBase<const_iterator,size_type,difference_type>
    {
        typedef IteratorBase<const_iterator,size_type,difference_type> BaseT;
    public:
        const_iterator() : mpHBV(0) {}
        const_iterator( HugeBVec const* pHBV, size_type pos )
        : BaseT(pos), mpHBV(pHBV) {}

        // compiler-supplied copying and destructor are OK

        value_type operator*() const { return (*mpHBV)[pos()]; }
        value_type operator[]( difference_type idx ) const
        { return (*mpHBV)[pos()+idx]; }

    private:
        HugeBVec const* mpHBV;
    };

    class const_rc_iterator
    : public ItrTagBase,
      public IteratorBase<const_rc_iterator,size_type,difference_type>
    {
        typedef IteratorBase<const_rc_iterator,size_type,difference_type> BaseT;
    public:
        const_rc_iterator() : mpHBV(0), mLast(~0ul) {}
        const_rc_iterator( HugeBVec const* pHBV, size_type end, size_type pos )
        : BaseT(pos), mpHBV(pHBV), mLast(end-1) {}

        // compiler-supplied copying and destructor are OK

        value_type operator*() const
        { return (*mpHBV)[mLast-pos()] ^ 3; }

        value_type operator[]( difference_type idx ) const
        { return (*mpHBV)[mLast-(pos()+idx)] ^ 3; }

    protected:
        HugeBVec const* mpHBV;
        size_type mLast;
    };

    HugeBVec( char const* fileName );
    ~HugeBVec();

    size_type size() const { return mSize; }

    value_type operator[]( size_t idx ) const
    { AssertLt(idx,mSize); return (mpBuf[idx>>2] >> 2*(3-(idx&3)))&3; }

    const_iterator begin( size_type idx = 0 ) const
    { return const_iterator(this,idx); }
    const_iterator end() const
    { return const_iterator(this,size()); }
    const_rc_iterator rcbegin( size_type end, size_type idx = 0 ) const
    { return const_rc_iterator(this,end,idx); }
    const_rc_iterator rcend( size_type end ) const
    { return const_rc_iterator(this,end,end); }

    friend bool operator==( HugeBVec const& v1, HugeBVec const& v2 )
    { return v1.size()==v2.size() &&
             !memcmp(v1.mpBuf,v2.mpBuf,v1.physSize()); }

    friend bool operator!=( HugeBVec const& v1, HugeBVec const& v2 )
    { return !(v1 == v2); }

private:
    HugeBVec( HugeBVec const& ); // unimplemented -- no copying
    HugeBVec& operator=( HugeBVec const& ); // unimplemented -- no copying

    size_t physSize() const { return (mSize+3)/4; }

    unsigned char const* mpBuf;
    size_type mSize;
};

#endif /* FEUDAL_HUGEBVEC_H_ */
