///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file SerfVec.h
 * \author tsharpe
 * \date Jul 27, 2009
 *
 * \brief
 */
#ifndef FEUDAL_SERFVEC_H_
#define FEUDAL_SERFVEC_H_

#include "feudal/BinaryStreamTraits.h"
#include "feudal/SmallVec.h"
#include "feudal/Mempool.h"
#include "system/Assert.h"

template <class T>
class SerfVec : public SmallVec<T,MempoolAllocator<T> >
{
    typedef MempoolAllocator<T> Alloc;
    typedef SmallVec<T,MempoolAllocator<T> > BaseT;

public:
    typedef Alloc alloc_type;
    typedef typename BaseT::size_type size_type;

    SerfVec() {}

    explicit SerfVec( Alloc const& alloc ) : BaseT(alloc) {}

    explicit SerfVec( size_type n )
    : BaseT(n)
    {}

    explicit SerfVec( size_type n,
                      T const& exemplar,
                      size_type capacity = 0,
                      Alloc const& alloc = Alloc() )
    : BaseT(n,exemplar,capacity,alloc)
    {}

    template <class Itr>
    SerfVec( Itr first, Itr const& last,
             size_type capacity = 0,
             Alloc const& alloc = Alloc() )
    : BaseT(first,last,capacity,alloc) {}

    SerfVec( SerfVec const& that )
    : BaseT(that) {}

    // compiler-supplied destructor and copying is OK

    template <class A1>
    SerfVec& operator=( SmallVec<T,A1> const& that )
    { BaseT::operator=(that); return *this; }

    SerfVec& operator=( SerfVec const& that )
    { BaseT::operator=(that); return *this; }

    /// This is deprecated.  Ownership is handled by the allocator.
    bool SelfOwned() const { return false; }

    /// Deprecated:  Use clear().shrink_to_fit().
    void Reinitialize() { destroy(); }
    /// Deprecated:  Use clear().shrink_to_fit().
    void Blank() { destroy(); }
    /// Deprecated:  Use clear().shrink_to_fit().
    void destroy() { BaseT::clear(); BaseT::shrink_to_fit(); }

    /// Deprecated:  Use swap().
    void Swap( SerfVec& that ) { BaseT::swap(that); }

    /// Deprecated:  Use assign().
    SerfVec& SetToSubOf( SerfVec const& that, size_type pos, size_type len )
    { AssertLe(pos,that.size());
      AssertLe(len,that.size()-pos);
      if ( this != &that )
      { assign(that.begin(pos),that.begin(pos+len)); }
      else
      { erase(BaseT::begin(),BaseT::begin(pos));
        BaseT::resize(len); }
      return *this; }

    /// Deprecated:  Use an iterator-pair constructor with reverse iterators
    /// from that other SerfVec.
    SerfVec& SetToReverseOf( SerfVec const& that )
    { if ( this == &that ) ReverseMe();
      else assign(that.rbegin(),that.rend());
      return *this; }

    /// Deprecated:  Use std::reverse().
    SerfVec& ReverseMe()
    { using std::reverse;
      reverse(BaseT::begin(),BaseT::end());
      return *this; }

    // Deprecated: You don't need this anymore.
    size_type SizeOfDynamicData() const { return BaseT::size(); }

    friend SerfVec Cat( SerfVec const& left, SerfVec const& right )
    { SerfVec result(left.begin(),left.end(),left.size()+right.size());
      result.append(right.begin(),right.end());
      return result; }
};

template <class T>
struct Serializability< SerfVec<T> > : public SelfSerializable {};

template <class T>
void swap( SerfVec<T>& v1, SerfVec<T>& v2 )
{ v1.swap(v2); }

/// Deprecated:  Use std::reverse().
template <class T>
inline void ReverseThis( SerfVec<T>& v )
{ using std::reverse; reverse(v.begin(),v.end()); }

/// Deprecated:  Use iterator-pair constructor with reverse iterators.
template <class T>
inline SerfVec<T> Reverse( SerfVec<T> const& v )
{ return SerfVec<T>(v.rbegin(),v.rend()); }

/// Deprecated:  Use clear().shrink_to_fit().
template <class T>
inline void Destroy( SerfVec<T>& v )
{ v.destroy(); }

#endif /* FEUDAL_SERFVEC_H_ */
