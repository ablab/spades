///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file ID.h
 * \author tsharpe
 * \date Jun 15, 2010
 *
 * \brief An N-byte value to use as an ID.
 */
#ifndef SYSTEM_ID_H_
#define SYSTEM_ID_H_

#include "system/Assert.h"
#include "system/StaticAssert.h"
#include <cstddef>
#include <cstring>
#include <ostream>

/// An N-byte integer (N <= 8) that can be used as an ID.
/// Has a 1-byte alignment requirement.  (I.e., packs nicely.)
/// The largest value is reserved as a null ID.
template <unsigned N>
class ID
{
public:
    ID() { memset(mVal,-1,N); }

    explicit ID( size_t val ) { setVal(val); }

    // compiler-generated copying and destructor are OK

    size_t val() const
    { size_t result; memcpy(&result,mVal,N); return result & MASK; }

    void setVal( size_t val )
    { AssertNot(val & ~MASK); memcpy(mVal,&val,N); }

    bool isNull() const { return *this == gNull; }

    friend bool operator==( ID const& id1, ID const& id2 )
    { return !memcmp(id1.mVal,id2.mVal,N); }

    friend bool operator!=( ID const& id1, ID const& id2 )
    { return memcmp(id1.mVal,id2.mVal,N); }

    friend bool operator<( ID const& id1, ID const& id2 )
    { return id1.val() < id2.val(); }

    friend bool operator<=( ID const& id1, ID const& id2 )
    { return id1.val() <= id2.val(); }

    friend bool operator>( ID const& id1, ID const& id2 )
    { return id1.val() > id2.val(); }

    friend bool operator>=( ID const& id1, ID const& id2 )
    { return id1.val() >= id2.val(); }

    friend int compare( ID const& id1, ID const& id2 )
    { return compare(id1.val(),id2.val()); }

    friend std::ostream& operator<<( std::ostream& os, ID id )
    { return os << id.val(); }

private:
    ID( bool )
    { STATIC_ASSERT(N <= 8); assertLittleEndian(); memset(mVal,-1,N); }

    void assertLittleEndian()
    { long foo = 1;
      ForceAssertEq(*reinterpret_cast<char*>(&foo),1); }

    unsigned char mVal[N];
    static ID gNull;
    static size_t const MASK = (1ul << (8*N)) - 1ul;
};

template <unsigned N> ID<N> ID<N>::gNull(true);

#endif /* SYSTEM_ID_H_ */
