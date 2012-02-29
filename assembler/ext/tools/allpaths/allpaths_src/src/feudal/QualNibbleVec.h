/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file QualNibbleVec.h
 * \author iainm
 * \date Sep 1, 2009
 *
 * \brief Quality scores compressed to a half-byte.
 */
#ifndef QUAL_NIBBLE_VECTOR_H
#define QUAL_NIBBLE_VECTOR_H

#include "String.h"
#include "Qualvector.h"
#include "feudal/FieldVec.h"
#include "feudal/MasterVec.h"
#include "feudal/Mempool.h"
#include "system/Assert.h"
#include <algorithm>
#include <cstddef>
#include <iterator>

class QualNibbleVec : private FieldVec<4, MempoolAllocator<unsigned char> >
{
public:
  typedef allocator_type Alloc;
  typedef unsigned char value_type;
  typedef unsigned size_type;
  typedef FieldVec<4, MempoolAllocator<unsigned char> > BaseT;
  typedef std::ptrdiff_t difference_type;
  typedef std::iterator<std::random_access_iterator_tag,
                        value_type,
                        difference_type,
                        void,
                        value_type> ItrBase;
  struct QualMapper
  {
      value_type operator()( value_type val ) const
      { return deflate(val); }

      static value_type deflate( value_type val )
      { using std::min; return min(val/3,15); }
      static value_type inflate( value_type val )
      { return 3*val + 1; }
  };

  class iterator
  : public ItrBase,
    public IteratorBase<iterator,size_type,difference_type>
  {
  public:
      iterator() : mpContainer(0) {}
      iterator( QualNibbleVec* pContainer, size_type pos )
      : IteratorBase<iterator,size_type,difference_type>(pos),
        mpContainer(pContainer) {}

      // compiler-supplied copying and destructor are OK

      value_type operator*() const { return (*mpContainer)[this->pos()]; }

      value_type operator[]( difference_type idx ) const
      { return (*mpContainer)[this->pos()+idx]; }

      void set( value_type val )
      { mpContainer->set(this->pos(),val); }

  protected:
      QualNibbleVec* mpContainer;
  };

  class const_iterator
  : public ItrBase,
    public IteratorBase<const_iterator,size_type,difference_type>
  {
  public:
      const_iterator() : mpContainer(0) {}
      const_iterator( QualNibbleVec const* pContainer, size_type pos )
      : IteratorBase<const_iterator,size_type,difference_type>(pos),
        mpContainer(pContainer) {}

      // compiler-supplied copying and destructor are OK

      value_type operator*() const { return (*mpContainer)[this->pos()]; }
      value_type operator[]( difference_type idx ) const
      { return (*mpContainer)[this->pos()+idx]; }

  private:
      QualNibbleVec const* mpContainer;
  };

  class reverse_iterator
  : public ItrBase,
    public IteratorBase<reverse_iterator,size_type,difference_type>
  {
  public:
      reverse_iterator() : mpContainer(0), mLast(~0) {}
      reverse_iterator( QualNibbleVec* pContainer, size_type pos )
      : IteratorBase<reverse_iterator,size_type,difference_type>(pos),
        mpContainer(pContainer), mLast(pContainer->size()-1) {}

      // compiler-supplied copying and destructor are OK

      value_type operator*() const { return (*mpContainer)[mLast-this->pos()]; }

      value_type operator[]( difference_type idx ) const
      { return (*mpContainer)[mLast-(this->pos()+idx)]; }

      void set( value_type val )
      { mpContainer->set(mLast-this->pos(),val); }

  protected:
      QualNibbleVec* mpContainer;
      size_type mLast;
  };

  class const_reverse_iterator
  : public ItrBase,
    public IteratorBase<const_reverse_iterator,size_type,difference_type>
  {
  public:
      const_reverse_iterator() : mpContainer(0), mLast(~0) {}
      const_reverse_iterator( QualNibbleVec const* pContainer, size_type pos )
      : IteratorBase<const_reverse_iterator,size_type,difference_type>(pos),
        mpContainer(pContainer), mLast(pContainer->size()-1) {}

      // compiler-supplied copying and destructor are OK

      value_type operator*() const { return (*mpContainer)[mLast-this->pos()]; }

      value_type operator[]( difference_type idx ) const
      { return (*mpContainer)[mLast-(this->pos()+idx)]; }

  protected:
      QualNibbleVec const* mpContainer;
      size_type mLast;
  };

  //
  // Constructors
  //
  QualNibbleVec() {}
  QualNibbleVec( Alloc const& alloc ) : BaseT(alloc) {}

  // SetToSubOf constructor.
  QualNibbleVec( const QualNibbleVec& q, const size_type start, const size_type len )
  { AssertLe(start,q.size()); AssertLe(len,q.size()-start);
    assign(q.begin(start),q.begin(start+len)); }

  // sized constructor
  explicit QualNibbleVec( size_type sz, size_type extra = 0 )
  : BaseT(sz, value_type(), sz+extra) {}

  // Copy constructor
  QualNibbleVec( QualNibbleVec const& q )
  : BaseT(q) {}

  // Construct from an 8 bit qualvector
  explicit QualNibbleVec( qualvector const & qv )
  { assign(qv.begin(),qv.end(),QualMapper()); }

  QualNibbleVec& operator=( QualNibbleVec const& q )
  { BaseT::operator=(q); return *this; }

  //
  // iterators
  //

  iterator begin() { return iterator(this,0); }
  iterator begin( size_type idx )
  { AssertLe(idx,size()); return iterator(this,idx); }
  iterator end() { return iterator(this,size()); }

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

  size_type size() const { return BaseT::size(); }
  size_type allocSize() const { return BaseT::allocSize(); }
  void swap( QualNibbleVec& q ) { BaseT::swap(q); }
  QualNibbleVec& reserve( size_type nnn ) { BaseT::reserve(nnn); return *this; }
  QualNibbleVec& resize( size_type nnn, value_type vvv = 0 )
  { BaseT::resize(nnn,QualMapper::deflate(vvv)); return *this; }
  QualNibbleVec& clear() { BaseT::clear(); return *this; }
  QualNibbleVec& push_back( value_type vvv )
  { BaseT::push_back(QualMapper::deflate(vvv)); return *this; }
  QualNibbleVec& reverse() { BaseT::reverse(); return *this; }

  size_t writeFeudal( BinaryWriter& writer, void const** ppFixed ) const
  { return BaseT::writeFeudal(writer,ppFixed); }
  void readFeudal( BinaryReader& rdr, size_t varDataLen, void* pFixed )
  { BaseT::readFeudal(rdr,varDataLen,pFixed); }
  size_t writeBinary( BinaryWriter& writer ) const
  { return BaseT::writeBinary(writer); }
  void readBinary( BinaryReader& reader ) { BaseT::readBinary(reader); }
  static size_t externalSizeof() { return 0; }
  //
  // Converstions to and from 8 bit qualvector
  //

  // Initialize QualNibbleVec from an 8 bit qualvector
  void SetFromQualvector( const qualvector& q )
  { assign(q.begin(),q.end(),QualMapper()); }

  // Convert QualNibbleVec to an 8 bit qualvector
  qualvector GetQualvector() const
  { qualvector qv; qv.reserve(size());
    for (size_type i = 0; i < size(); ++i) qv.push_back( (*this)[i] );
    return qv; }


  //
  // Accessors
  //

  unsigned char operator[]( size_type idx ) const
  { return QualMapper::inflate(BaseT::operator[](idx)); }

  unsigned char getRaw( size_type idx ) const
  { return BaseT::operator [](idx); }

  void set( size_type idx, unsigned char value)
  { BaseT::set(idx, QualMapper::deflate(value) ); }

  void setRaw( size_type idx, unsigned char value )
  { BaseT::set(idx,value); }

  void ReverseMe() { reverse(); }

  /// Replaces each quality score with the minimum quality score in the range
  /// idx-radius to idx+radius.
  QualNibbleVec& squash( unsigned radius );
};

SELF_SERIALIZABLE(QualNibbleVec);

typedef MasterVec<QualNibbleVec> VecQualNibbleVec;
typedef VecQualNibbleVec QualNibbleVecVec;
typedef VecQualNibbleVec vecqnibvec;

typedef QualNibbleVec qnibble ;

/// Store a VecQualNibbleVec in qualb format.
void WriteAll(QualNibbleVecVec const& quals, String const& fn);

/// Load a VecQualNibbleVec from a qualb file.
void LoadQualNibbleVec( const String & fn, VecQualNibbleVec * quals );

#endif
