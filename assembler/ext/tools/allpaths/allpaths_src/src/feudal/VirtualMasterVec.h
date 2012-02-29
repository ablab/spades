///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file VirtualMasterVec.h
 * \author tsharpe
 * \date Aug 18, 2009
 *
 * \brief A class that gives you sequential, read-only access to a feudal file.
 *
 * The memory footprint is very small, and sequential access is very efficient,
 * involving minimal I/O.
 *
 * There is a random access capability as well, but it requires a somewhat
 * expensive file seek and buffer purge and reload.  Nonetheless, the random
 * access capability might be quite suitable for sparse access.
 *
 * The template class T must be a type that was written to a feudal file, so it
 * will necessarily have the right interface, namely a copy-constructor and a
 * constructor that takes void*'s to the raw data:
 * T( void* varData, ulong varDataNBytes, void* fixedData, A const& allocator );
 */
#ifndef FEUDAL_VIRTUALMASTERVEC_H_
#define FEUDAL_VIRTUALMASTERVEC_H_

#include "feudal/FeudalFileReader.h"
#include "feudal/Oob.h"
#include <cstddef>

template<class T>
class VirtualMasterVec
{
public:
    typedef T const value_type;
    typedef unsigned long size_type;

    class Itr
    : public std::iterator<std::input_iterator_tag,T,std::ptrdiff_t>
    {
    public:
        Itr( VirtualMasterVec<T> const& vec, size_type idx )
        : mVec(vec), mIdx(idx), mLoadedIdx(~0ul) {}

        Itr( Itr const& that )
        : mVec(that.mVec), mIdx(that.mIdx), mLoadedIdx(~0UL) {}

        Itr& operator=( Itr const& that ) { mIdx = that.mIdx; return *this; }

        // compiler-supplied destructor is OK

        T const& operator*() const
        { if ( mLoadedIdx != mIdx ) { mRef = mVec[mIdx]; mLoadedIdx = mIdx; }
          return mRef; }
        T const* operator->() const { return &operator*(); }
        bool operator==( Itr const& that ) const { return mIdx == that.mIdx; }
        bool operator!=( Itr const& that ) const { return mIdx != that.mIdx; }
        Itr& operator++() { mIdx += 1; return *this; }
        Itr operator++(int) { Itr tmp(*this); mIdx += 1; return tmp; }

        bool operator<( Itr const& itr ) { return mIdx < itr.mIdx; }
        std::ptrdiff_t operator-( Itr const& itr ) { return mIdx - itr.mIdx; }
        Itr operator+( std::ptrdiff_t inc )
        { Itr tmp(*this); tmp.mIdx += inc; return tmp; }

    private:
        VirtualMasterVec const& mVec;
        size_type mIdx;
        mutable T mRef;
        mutable size_type mLoadedIdx;
    };

    typedef Itr const_iterator;

    VirtualMasterVec( char const* path ) : mFFR(path) {}

    /// Construct from something with a c_str() member (like string or String)
    template <class C>
    explicit VirtualMasterVec( C const& path,
                                    char const*(C::*)() const=&C::c_str )
    : mFFR(path.c_str()) {}

    VirtualMasterVec( VirtualMasterVec<T> const& that ) : mFFR(that.mFFR) {}
    // compiler-supplied destructor is OK
    VirtualMasterVec& operator=( VirtualMasterVec<T> const& that )
    { if ( this != &that ) { mFFR = that.mFFR; } return *this; }

    Itr begin() const { return Itr(*this,0); }
    Itr begin( size_type idx ) const { return Itr(*this,idx); }
    Itr end() const { return Itr(*this,mFFR.getNElements()); }

    T const front() const { return obj(0); }
    T const back() const { return obj(size()-1); }
    T const operator[]( size_type idx ) const { return obj(idx); }
    T const at( size_type idx ) const
    { if ( idx >= size() )
      { OutOfBoundsReporter::oob("VirtualMasterVec",idx,size()); }
      return obj(idx); }

    size_type size() const { return mFFR.getNElements(); }

    /// total number of T::value_type's across all elements.
    /// this will blow up if T::value_type doesn't have a fixed external size
    size_t sizeSum() const
    { typename T::value_type* ppp = static_cast<typename T::value_type*>(0);
      size_t eleSiz = BinaryReader::externalSizeof(ppp);
      AssertNe(eleSiz,0u); return mFFR.getDataLenTotal()/eleSiz; }

    bool empty() const { return size() == 0; }

    size_t getMappedLen() const { return mFFR.getMappedLen(); }

private:
    T const obj( size_type idx ) const
    { AssertLt(idx,size());
      T result;
      result.readFeudal(mFFR.getData(idx),mFFR.getDataLen(idx),
                         mFFR.getFixedData(idx,T::fixedDataLen()));
      return result; }

    mutable FeudalFileReader mFFR;
};

#endif // FEUDAL_VIRTUALMASTERVEC_H_
