///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Iterator.h
 * \author tsharpe
 * \date Jun 9, 2009
 *
 * \brief Random-access iterators for vector-like classes.
 */
#ifndef FEUDAL_ITERATOR_H_
#define FEUDAL_ITERATOR_H_

#include <cstddef>
#include <iterator>

/// An abstract iterator position.
template <class S>
class IteratorPosition
{
public:
    explicit IteratorPosition( S pos ) : mPos(pos) {}

    // compiler-supplied copying and destructor are OK

    S pos() const { return mPos; }

    friend bool operator==(IteratorPosition const&i1,IteratorPosition const&i2)
    { return i1.mPos == i2.mPos; }

    friend bool operator!=(IteratorPosition const&i1,IteratorPosition const&i2)
    { return i1.mPos != i2.mPos; }

    friend bool operator<(IteratorPosition const&i1,IteratorPosition const&i2)
    { return i1.mPos < i2.mPos; }

    friend bool operator<=(IteratorPosition const&i1,IteratorPosition const&i2)
    { return i1.mPos <= i2.mPos; }

    friend bool operator>(IteratorPosition const&i1,IteratorPosition const&i2)
    { return i1.mPos > i2.mPos; }

    friend bool operator>=(IteratorPosition const&i1,IteratorPosition const&i2)
    { return i1.mPos >= i2.mPos; }

protected:
    S mPos;
};

/// A piece of a bidirectional iterator for vector-like classes.
/// This piece implements position increments and decrements.
template <class X,class S,class D>
class IteratorBiDiBase : public IteratorPosition<S>
{
public:
    IteratorBiDiBase() : IteratorPosition<S>(0) {}
    explicit IteratorBiDiBase( S pos ) : IteratorPosition<S>(pos) {}

    // compiler-supplied copying and destructor are OK

    X& operator++() { thisX().fwd(1); return thisX(); }
    X operator++(int) { X tmp(thisX()); thisX().fwd(1); return tmp; }
    X& operator--() { thisX().bwd(1); return thisX(); }
    X operator--(int) { X tmp(thisX()); thisX().bwd(1); return tmp; }

    void fwd( D diff ) { this->mPos += diff; }
    void bwd( D diff ) { this->mPos -= diff; }

private:
    X& thisX() { return *static_cast<X*>(this); }
    X const& thisX() const { return *static_cast<X const*>(this); }
};


/// A piece of a random-access iterator for vector-like classes.
/// This piece implements random position movement.
template <class X,class S,class D>
class IteratorBase : public IteratorBiDiBase<X,S,D>
{
    typedef IteratorBiDiBase<X,S,D> Base;
public:
    IteratorBase() : Base(0) {}
    explicit IteratorBase( S pos ) : Base(pos) {}

    // compiler-supplied copying and destructor are OK

    X& operator+=( D diff ) { thisX().fwd(diff); return thisX(); }

    X operator+( D diff ) const
    { X tmp(thisX()); tmp.thisX().fwd(diff); return tmp; }

    X& operator-=( D diff ) { thisX().bwd(diff); return thisX(); }

    X operator-( D diff ) const
    { X tmp(thisX()); tmp.thisX().bwd(diff); return tmp; }

    D operator-( IteratorPosition<S> const& that ) const
    { return this->mPos - that.pos(); }

    void fwd( D diff ) { this->mPos += diff; }
    void bwd( D diff ) { this->mPos -= diff; }

private:
    X& thisX() { return *static_cast<X*>(this); }
    X const& thisX() const { return *static_cast<X const*>(this); }
};

/*
 * Iterator position differences.
 */

template <class X,class S,class D>
X operator+( D diff, IteratorBase<X,S,D> const& itr )
{ return itr+diff; }

template <class T>
class FwdIterator
: public std::iterator<std::random_access_iterator_tag,T>,
  public IteratorBase<FwdIterator<T>,T*,std::ptrdiff_t>
{
    typedef IteratorBase<FwdIterator<T>,T*,std::ptrdiff_t> Base;
public:
    FwdIterator() {}
    FwdIterator( T* pos ) : Base(pos) {}

    // compiler-supplied copying, destruction are OK

    T& operator*() const { return *this->mPos; }
    T* operator->() const { return this->mPos; }
    T& operator[]( std::ptrdiff_t diff ) const { return this->mPos[diff]; }
};

template <class T>
class FwdConstIterator
: public std::iterator<std::random_access_iterator_tag,T,std::ptrdiff_t,
                          T const*,T const&>,
  public IteratorBase<FwdConstIterator<T>,T const*,std::ptrdiff_t>
{
    typedef IteratorBase<FwdConstIterator<T>,T const*,std::ptrdiff_t> Base;
public:
    FwdConstIterator() {}

    FwdConstIterator( T const* pos ) : Base(pos) {}

    FwdConstIterator( IteratorPosition<T const*> const& pos )
    : Base(pos.pos()) {}

    FwdConstIterator( IteratorPosition<T*> const& pos ) : Base(pos.pos()) {}

    // compiler-supplied copying, destruction are OK

    T const& operator*() const { return *this->mPos; }
    T const* operator->() const { return this->mPos; }
    T const& operator[]( std::ptrdiff_t diff ) const
    { return this->mPos[diff]; }
};

template <class T>
class RevIterator
: public std::iterator<std::random_access_iterator_tag,T>,
  public IteratorBase<RevIterator<T>,unsigned long,std::ptrdiff_t>
{
    typedef IteratorBase<RevIterator<T>,unsigned long,std::ptrdiff_t> Base;
public:
    RevIterator() {}

    RevIterator( unsigned long idx, T* pEnd ) : Base(idx), mpEnd(pEnd) {}

    // compiler-supplied copying, destruction are OK

    T& operator*() const { return mpEnd[-this->mPos-1]; }
    T* operator->() const { return mpEnd-this->mPos-1; }
    T& operator[]( std::ptrdiff_t diff ) const
    { return mpEnd[-this->mPos-diff-1]; }

private:
    T* mpEnd;
};

template <class T>
class RevConstIterator
: public std::iterator<std::random_access_iterator_tag,T,std::ptrdiff_t,
                          T const*,T const&>,
  public IteratorBase<RevConstIterator<T>,unsigned long,std::ptrdiff_t>
{
    typedef IteratorBase<RevConstIterator<T>,unsigned long,std::ptrdiff_t> Base;
public:
    RevConstIterator() {}

    RevConstIterator( unsigned long idx, T const* pEnd )
    : Base(idx), mpEnd(pEnd) {}

    RevConstIterator( IteratorPosition<unsigned long> const& pos )
    : Base(pos.pos()) {}

    // compiler-supplied copying, destruction are OK

    T const& operator*() const { return mpEnd[-this->mPos-1]; }
    T const* operator->() const { return mpEnd-this->mPos-1; }

    T const& operator[]( std::ptrdiff_t diff ) const
    { return mpEnd[-this->mPos-diff-1]; }

private:
    T const* mpEnd;
};

#endif /* FEUDAL_ITERATOR_H_ */
