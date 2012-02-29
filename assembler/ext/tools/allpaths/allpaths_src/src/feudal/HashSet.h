///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file HashSet.h
 * \author tsharpe
 * \date Aug 20, 2010
 *
 * \brief Big, low-memory-overhead hash table.
 */
#ifndef FEUDAL_HASHSET_H
#define FEUDAL_HASHSET_H

#include "feudal/BinaryStream.h"
#include "feudal/Iterator.h"
#include "math/PowerOf2.h"
#include "system/Assert.h"
#include "system/SpinLockedData.h"
#include <algorithm>
#include <cstddef>
#include <cstring>

/// Hopscotch hash set needs an abstract way of producing an array of values
/// that can be used when T doesn't have a default constructor (e.g., Kmers).
/// This implementation is for the simple case where you can use new and delete.
template <class T>
class TFactory
{
public:
    T* create( size_t nTs ) const { return new T[nTs]; }
    void destroy( T* pTs, size_t nTs ) const { delete [] pTs; }
};

// Ignore this little helper class.
// It merely subclasses some likely-to-be-empty classes so they don't take up
// extra space in the hash set classes.
/// Template args as for HopscotchHashSet, which follows.
template <class T, class H, class C = std::equal_to<T>, class F=TFactory<T> >
class HCFLock : public H, public C, public F, public SpinLockedData
{
public:
    typedef typename H::argument_type const& key_type;

    explicit HCFLock( H const& h= H(), C const& c = C(), F const& f = F() )
    : H(h), C(c), F(f) {}

    HCFLock( HCFLock const& that )
    : H(that), C(that), F(that), SpinLockedData() {}

    // compiler-supplied destructor is OK

    size_t hash( key_type val ) const
    { return static_cast<H const&>(*this)(val); }
    bool compare( key_type v1, key_type v2 ) const
    { return static_cast<C const&>(*this)(v1,v2); }

    friend void swap( HCFLock& hcf1, HCFLock& hcf2 )
    { using std::swap;
      swap(static_cast<H&>(hcf1),static_cast<H&>(hcf2));
      swap(static_cast<C&>(hcf1),static_cast<C&>(hcf2));
      swap(static_cast<F&>(hcf1),static_cast<F&>(hcf2)); }
      // lock not swapped, and that's on purpose

private:
    HCFLock& operator=( HCFLock const& ); // unimplemented -- no copying
};

class BktStatus
{
public:
    enum Status { EMPTY, HEAD, SQUATTER };

    // compiler-supplied default constructor, copying and destructor are OK

    Status getStatus() const { return static_cast<Status>(mStatus); }
    void setStatus( Status status ) { mStatus = status; }

    unsigned getOffset() const { return mOffset; }
    void setOffset( unsigned offset )
    { AssertLe(offset,maxOffset()); mOffset = offset; }

    bool hasNext() const { return mOffset; }
    BktStatus const* next() const { return this-mOffset; }
    BktStatus* next() { return this-mOffset; }

    static unsigned maxOffset() { return (1ul << 6) - 1ul; }

private:
    unsigned char mStatus : 2;
    unsigned char mOffset : 6;
};

TRIVIALLY_SERIALIZABLE(BktStatus);

// fwd decl for friend
template <class T, class H, class C, class F> class HashSet;

/// A fixed-size hash set.
/// The hopscotch algorithm provides high load factors without getting in a jam.
/// This implementation is optimized for low memory overhead: just a single
/// byte per entry.
///
/// T, the value_type, must be assignable.  It must also be default
/// constructable, unless you supply a custom factory.
///
/// H is a hashing functor.
/// It transforms a T into a size_t, it's hash value:
/// size_t H( T const& ) const;
/// It must have a typedef called argument_type.
///
/// C is a comparison functor.
/// It returns true if two T's are equal:
/// bool C( T const&, T const& ) const;
///
/// F is a factory for an array of Ts.
///
/// This class is NOT thread-safe.  However, it has a lock and unlock method,
/// so you can single-thread it, if you want.  That's what HashSet does.
template <class T, class H, class C = std::equal_to<T>, class F=TFactory<T> >
class HopscotchHashSet
{
public:
    typedef T value_type;
    typedef value_type const* const_pointer;
    typedef value_type const& const_reference;
    typedef unsigned size_type;
    typedef std::ptrdiff_t difference_type;
    typedef H hash_func_type;
    typedef C comp_func_type;
    typedef F factory_type;
    typedef typename H::argument_type const& key_type;
    friend class HashSet<T,H,C,F>;

    class NoRoomException : public std::exception
    {
    public:
        NoRoomException( char const* msg ) : mMsg(msg) {}
        char const* what() const throw() { return mMsg; }
    private:
        char const* mMsg;
    };

    // this is a const_iterator over the value_type entries
    class const_iterator
    : public std::iterator<std::bidirectional_iterator_tag,value_type const>,
      public IteratorBiDiBase<const_iterator,const_pointer,difference_type>
    {
        typedef IteratorBiDiBase<const_iterator,const_pointer,difference_type> Base;
    public:
        const_iterator() : Base(0), mpHHS(0) {}

        const_iterator( const_pointer pCur, HopscotchHashSet const* pHS )
        : Base(pCur-1), mpHHS(pHS)
        { fwd(1); }

        // compiler-supplied copying and destructor are OK

        const_reference operator*() const { return *Base::pos(); }
        const_pointer operator->() const { return Base::pos(); }

        void fwd( difference_type diff )
        { if ( diff < 0 ) bwd(-diff);
          else
            while ( diff-- )
              do { Base::fwd(1); }
              while ( mpHHS->isEmptyPosition(Base::pos()) ); }

        void bwd( difference_type diff )
        { if ( diff < 0 ) fwd(-diff);
          else
            while ( diff-- )
              do { Base::bwd(1); }
              while ( mpHHS->isEmptyPosition(Base::pos()) ); }

    private:
        HopscotchHashSet const* mpHHS;
    };

    HopscotchHashSet( size_type cap, HCFLock<T,H,C,F>& hcf )
    : mCapacity(adjustCapacity(cap)), mSize(0),
      mBuckets(hcf.create(capacity())), mBktInfo(new BktStatus[capacity()]),
      mHCF(hcf)
    { memset(mBktInfo,0,capacity()); }

    HopscotchHashSet( HopscotchHashSet const& that )
    : mCapacity(that.capacity()), mSize(0),
      mBuckets(that.mHCF.create(capacity())),
      mBktInfo(new BktStatus[capacity()]), mHCF(that.mHCF)
    { memset(mBktInfo,0,capacity());
      *this = that; }

    ~HopscotchHashSet()
    { mHCF.destroy(mBuckets,capacity()); delete [] mBktInfo; }

    HopscotchHashSet& operator=( HopscotchHashSet const& that )
    { if ( this != &that )
      { clear();
        BktStatus const* start = that.mBktInfo;
        BktStatus const* end = start + that.capacity();
        for ( BktStatus const* itr(start); itr != end; ++itr )
          if ( itr->getStatus() == BktStatus::HEAD )
            add(that.mBuckets[itr-start]);
        for ( BktStatus const* itr(start); itr != end; ++itr )
          if ( itr->getStatus() == BktStatus::SQUATTER )
            add(that.mBuckets[itr-start]); }
      return *this; }

    size_type capacity() const { return mCapacity; }
    size_type size() const { return mSize; }

    const_iterator begin() const { return const_iterator(mBuckets,this); }
    const_iterator end() const
    { return const_iterator(mBuckets+capacity(),this); }

    /// Returns true if the value is in the set.
    bool contains( key_type val ) const { return lookup(val); }

    /// Returns pointer to val, or null pointer if val is not in set.
    const_pointer lookup( key_type val ) const
    { return lookup(val,mHCF.hash(val)); }

    /// Returns true if the value was added (otherwise, it was already present).
    bool add( key_type val ) throw(NoRoomException)
    { return add(val,mHCF.hash(val)); }

    /// Returns a reference to the value.  Inserts value, if not present.
    const_reference operator[]( key_type val ) throw(NoRoomException)
    { size_type idx = mHCF.hash(val) % capacity();
      const_pointer pBkt = find(val,idx);
      if ( pBkt ) return *pBkt; // EARLY RETURN!
      return (findInsertSlot(val,idx) = val); }

    /// Removes the value from the set.  Returns false if value not present.
    bool remove( key_type val )
    { return remove(val,mHCF.hash(val)); }

    /// Removes all values from the set.
    HopscotchHashSet& clear()
    { memset(mBktInfo,0,capacity()); mSize = 0; return *this; }

    void lock() { mHCF.lock(); }
    void unlock() { mHCF.unlock(); }

    bool isEmptyPosition( const_pointer pBkt ) const
    { size_type idx = pBkt - mBuckets;
      return idx < capacity() && mBktInfo[idx].getStatus()==BktStatus::EMPTY; }

    size_t writeBinary( BinaryWriter& writer ) const
    { size_t len = writer.write(mSize);
      len += writer.write(mBktInfo,mBktInfo+mCapacity);
      return len+writer.write(mBuckets,mBuckets+mCapacity); }

    void readBinary( BinaryReader& reader )
    { reader.read(&mSize);
      reader.read(mBktInfo,mBktInfo+mCapacity);
      reader.read(mBuckets,mBuckets+mCapacity); }

    static size_t externalSizeof() { return 0; }

    /// Swaps set contents.
    friend void swap( HopscotchHashSet& set1, HopscotchHashSet& set2 )
    { using std::swap;
      swap(set1.mCapacity,set2.mCapacity);
      swap(set1.mSize,set2.mSize);
      swap(set1.mBuckets,set2.mBuckets);
      swap(set1.mBktInfo,set2.mBktInfo);
      swap(set1.mHCF,set2.mHCF); }
      // note:  lock is not swapped, and that's on purpose

protected:
    // these protected methods are used by the friendly HashSet, which follows

    const_pointer lookup( key_type val, size_t hash ) const
    { return find(val,hash%capacity()); }

    bool add( key_type val, size_t hash ) throw(NoRoomException)
    { size_type idx = hash%capacity();
      if ( find(val,idx) ) return false;
      findInsertSlot(idx) = val; return true; }

    bool remove( key_type val, size_t hash )
    { // Find the head slot, and make sure it's actually a chain head. (If it's
      // empty or an overflow slot, there are no keys with the right hash.)
      BktStatus* pInfo = mBktInfo + hash%capacity();
      if ( pInfo->getStatus() != BktStatus::HEAD ) return false;// EARLY RETURN!
      BktStatus* pInfoPrev = 0;
      value_type* pBkt;
      // Find the right slot, where the value compares equal to the one sought.
      while ( true )
      { pBkt = &mBuckets[pInfo - mBktInfo];
        if ( mHCF.compare(val,*pBkt) ) break;
        if ( !pInfo->hasNext() ) return false; // EARLY RETURN!
        pInfoPrev = pInfo;
        pInfo = next(pInfo); }

      // At this point, pInfo and pBkt point to the right places for the value
      // to be deleted, and pInfoPrev is either the previous entry in the chain,
      // or null (when the value was found at the chain head).
      if ( pInfo->hasNext() ) copyAndKillEndOfChain(pInfo,pBkt);
      else // else at end of chain
      { if ( pInfoPrev ) pInfoPrev->setOffset(0);
        pInfo->setStatus(BktStatus::EMPTY); }
      mSize -= 1;
      return true; }

    const_reference insertK( key_type val, size_t hash ) throw(NoRoomException)
    { return (findInsertSlot(hash%capacity()) = val); }

    void insertV( const_reference val, size_t hash ) throw(NoRoomException)
    { findInsertSlot(hash%capacity()) = val; }

private:
    BktStatus* next( BktStatus* pInfo )
    { pInfo = pInfo->next();
      if ( pInfo < mBktInfo ) pInfo += capacity();
      return pInfo; }

    BktStatus const* next( BktStatus const* pInfo ) const
    { pInfo = pInfo->next();
      if ( pInfo < mBktInfo ) pInfo += capacity();
      return pInfo; }

    size_type distance( BktStatus const* pInfo, BktStatus const* pInfo2 ) const
    { return pInfo2 < pInfo ? pInfo-pInfo2 : pInfo-pInfo2+capacity(); }

    // Walk chain to find value.
    const_pointer find( key_type val, size_type offset ) const
    { // Check home slot.
      const_pointer pBkt = 0;
      BktStatus const* pInfo = mBktInfo + offset;
      if ( pInfo->getStatus() == BktStatus::HEAD )
      { pBkt = &mBuckets[pInfo - mBktInfo];
        if ( !mHCF.compare(val,*pBkt) )
        { pBkt = 0;
          // Check the rest of the chain, if any.
          while ( pInfo->hasNext() )
          { pInfo = next(pInfo);
            const_reference bkt = mBuckets[pInfo - mBktInfo];
            if ( mHCF.compare(val,bkt) ) {pBkt = &bkt; break; } } } }
      return pBkt; }

    void copyAndKillEndOfChain( BktStatus* pInfo, value_type* pBkt )
    { Assert(pInfo->hasNext());
      BktStatus* pInfoPrev;
      do { pInfoPrev = pInfo; pInfo = next(pInfo); }
      while ( pInfo->hasNext() );
      *pBkt = mBuckets[pInfo-mBktInfo];
      pInfoPrev->setOffset(0);
      pInfo->setStatus(BktStatus::EMPTY); }

    // Find a place for a value known not to be present.  Idx is the chain head.
    value_type& findInsertSlot( size_type idx )
      throw(NoRoomException)
    { if ( size() == capacity() )
        throw NoRoomException("capacity exceeded");
      BktStatus* pInfo = mBktInfo + idx;
      switch ( pInfo->getStatus() )
      {case BktStatus::SQUATTER: // move squatter, then use now-empty slot
         bump(pInfo);
         // no break statement -- allow flow-through to next case
       case BktStatus::EMPTY: // just use the empty slot
         pInfo->setStatus(BktStatus::HEAD);
         break;
       case BktStatus::HEAD: // extend chain
       { while ( pInfo->hasNext() ) pInfo = next(pInfo);
         BktStatus* pInfo2 = findEmpty(pInfo);
         pInfo2->setStatus(BktStatus::SQUATTER);
         pInfo->setOffset(distance(pInfo,pInfo2));
         pInfo = pInfo2;
         break; } }
      mSize += 1;
      return mBuckets[pInfo-mBktInfo]; }

    // move a SQUATTER key elsewhere
    void bump( BktStatus* pBump ) throw(NoRoomException)
    { AssertEq((int)pBump->getStatus(),(int)BktStatus::SQUATTER);

      // search for the entry that points to pBump
      BktStatus* pPrev = pBump;
      size_type off = 0;
      do
      { if ( ++pPrev >= mBktInfo+capacity() ) pPrev -= capacity(); }
      while ( ++off != pPrev->getOffset() );

      // loop for weird edge case where pEmpty gets picked too far from
      // its successor-to-be
      BktStatus* pEmpty;
      size_type distToEmpty;
      size_type distToSucc;
      while ( true ) // loop for weird edge case
      { pEmpty = findEmpty(pPrev);
        distToEmpty = distance(pPrev,pEmpty);
        distToSucc = pPrev->getOffset() + pBump->getOffset();
        if ( !pBump->hasNext() ||
             distToEmpty > pPrev->getOffset() ||
             distToSucc-distToEmpty <= BktStatus::maxOffset() ) break;
        // handle difficult case where pEmpty is too far from its successor
        copyAndKillEndOfChain(pBump,&mBuckets[pEmpty-mBktInfo]);
        pEmpty->setStatus(BktStatus::SQUATTER);
        // link pEmpty into chain between pPrev and pBump
        pEmpty->setOffset(pPrev->getOffset()-distToEmpty);
        pPrev->setOffset(distToEmpty);
        pPrev = pEmpty; }

      // move pBump's value to pEmpty
      move(pPrev,pBump,pEmpty,distToEmpty,distToSucc); }

    void move( BktStatus* pPrev, BktStatus* pBump, BktStatus* pEmpty,
               size_type distToEmpty, size_type distToSucc )
    { // move pBump's value to pEmpty
      mBuckets[pEmpty-mBktInfo] = mBuckets[pBump-mBktInfo];
      pEmpty->setStatus(BktStatus::SQUATTER);

      // fix up chain, as necessary
      if ( pBump->hasNext() )
      { if ( distToSucc > distToEmpty )
          pEmpty->setOffset(distToSucc-distToEmpty);
        else
        { pPrev->setOffset(distToSucc); // snip pBump out of chain
          // walk down chain to the end, or until pEmpty is between pPrev
          // and its successor
          do { distToEmpty -= pPrev->getOffset(); pPrev = next(pPrev); }
          while ( pPrev->hasNext() && pPrev->getOffset() < distToEmpty );
          if ( pPrev->hasNext() )
            pEmpty->setOffset(pPrev->getOffset()-distToEmpty); } }
      pPrev->setOffset(distToEmpty);

      // pBump is now free
      pBump->setStatus(BktStatus::EMPTY);
      pBump->setOffset(0); }

    // find empty slot within range of pInfo, or create one by hopscotching
    BktStatus* findEmpty( BktStatus* pInfo ) throw(NoRoomException)
    { AssertLt(mSize,capacity());
      BktStatus* pEmpty = pInfo;
      do { if ( --pEmpty < mBktInfo ) pEmpty += capacity(); }
      while ( pEmpty->getStatus() != BktStatus::EMPTY );
      while ( distance(pInfo,pEmpty) > BktStatus::maxOffset() )
          pEmpty = hopScotch(pEmpty);
      return pEmpty; }

    BktStatus* hopScotch( BktStatus* pEmpty ) throw (NoRoomException)
    { AssertEq((int)pEmpty->getStatus(),(int)BktStatus::EMPTY);
      BktStatus* pPrev = pEmpty + BktStatus::maxOffset();
      if ( pPrev >= mBktInfo + capacity() ) pPrev -= capacity();
      size_type distToEmpty = BktStatus::maxOffset();
      // for each entry that's in range of the free slot
      while ( pPrev != pEmpty )
      { // if pPrev's successor follows pEmpty (a move in the right direction)
        if ( pPrev->hasNext() && pPrev->getOffset() < distToEmpty )
        { BktStatus* pBump = next(pPrev);
          size_type distToSucc = pPrev->getOffset()+pBump->getOffset();
          // if the chain can be made whole after swapping pBump and pEmpty
          if ( !pBump->hasNext() || distToEmpty > pPrev->getOffset() ||
                  distToSucc-distToEmpty <= BktStatus::maxOffset() )
          { move(pPrev,pBump,pEmpty,distToEmpty,distToSucc);
            return pBump; } } // EARLY RETURN!
        distToEmpty -= 1;
        if ( --pPrev < mBktInfo ) pPrev += capacity(); }
      throw NoRoomException("hopscotching failed"); }

    static size_type adjustCapacity( size_type capacity )
    { if ( capacity < 2*BktStatus::maxOffset()+1 )
        capacity = 2*BktStatus::maxOffset()+1;
      capacity |= 1;
      while ( !capacity%3 || !capacity%5 || !capacity%7 || !capacity%11 )
        capacity += 2;
      return capacity; }

    size_type mCapacity;
    size_type mSize;
    value_type* mBuckets;
    BktStatus* mBktInfo;
    HCFLock<T,H,C,F> mHCF;
};

template <class T, class H, class C, class F>
struct Serializability< HopscotchHashSet<T,H,C,F> >
: public SelfSerializable {};

/// A gracefully-growable HashSet.
/// Splits fixed-size hopscotch hash sets as it grows.  This means that only a
/// fraction of the keys might need to be reorganized for any given call to
/// add() or operator[], which is smoother than reorganizing the whole
/// enchilada, as is typical in size-doubling schemes.  It also means that
/// memory use is in a few, large chunks, rather than requiring one giant
/// contiguous space.
/// Template args as for HopscotchHashSet.
///
/// Thread-safe for multiple writers.  Thread-safe for multiple readers.
/// NOT thread-safe for a mix of readers and writers.
template <class T, class H, class C=std::equal_to<T>, class F=TFactory<T> >
class HashSet
{
public:
    typedef T value_type;
    typedef value_type const* const_pointer;
    typedef value_type const& const_reference;
    typedef size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef H hash_func_type;
    typedef C comp_func_type;
    typedef F factory_type;
    typedef typename H::argument_type const& key_type;
    typedef HopscotchHashSet<T,H,C,F> HHS;
    typedef HHS*volatile* PPHHS;
    typedef typename HHS::NoRoomException NoRoomException;
    typedef typename HHS::const_iterator ICItr;

    // This is an iterator over the HopscotchHashSets.
    // Sorry to expose this bit of internal organization, but it seems simpler.
    // To iterate over entries, you'll need a nested for loop.
    class iterator
    : public std::iterator<std::bidirectional_iterator_tag,HHS>,
      public IteratorBiDiBase<iterator,PPHHS,difference_type>
    {
        typedef IteratorBiDiBase<iterator,PPHHS,difference_type> Base;
    public:
        iterator() : Base(0), mpHS(0) {}

        iterator( PPHHS pCur, HashSet* pHS )
        : Base(pCur-1), mpHS(pHS)
        { fwd(1); }

        // compiler-supplied copying and destructor are OK

        HHS& operator*() const { return **Base::pos(); }
        HHS* operator->() const { return *Base::pos(); }

        void fwd( difference_type diff )
        { if ( diff < 0 ) bwd(-diff);
          else
            while ( diff-- )
              do { Base::fwd(1); }
              while ( mpHS->isEmptyPosition(Base::pos()) ); }

        void bwd( difference_type diff )
        { if ( diff < 0 ) fwd(-diff);
          else
            while ( diff-- )
              do { Base::bwd(1); }
              while ( mpHS->isEmptyPosition(Base::pos()) ); }

    private:
        HashSet* mpHS;
    };

    // This is a const_iterator over the (const) HopscotchHashSets.
    // Sorry to expose this bit of internal organization, but it seems simpler.
    // To iterate over entries, you'll need a nested for loop.
    class const_iterator
    : public std::iterator<std::bidirectional_iterator_tag,HHS const>,
      public IteratorBiDiBase<const_iterator,PPHHS,difference_type>
    {
        typedef IteratorBiDiBase<const_iterator,PPHHS,difference_type> Base;
    public:
        const_iterator() : Base(0), mpHS(0) {}

        const_iterator( PPHHS pCur, HashSet const* pHS )
        : Base(pCur-1), mpHS(pHS)
        { fwd(1); }

        // compiler-supplied copying and destructor are OK

        HHS const& operator*() const { return **Base::pos(); }
        HHS const* operator->() const { return *Base::pos(); }

        void fwd( difference_type diff )
        { if ( diff < 0 ) bwd(-diff);
          else
            while ( diff-- )
              do { Base::fwd(1); }
              while ( mpHS->isEmptyPosition(Base::pos()) ); }

        void bwd( difference_type diff )
        { if ( diff < 0 ) fwd(-diff);
          else
            while ( diff-- )
              do { Base::bwd(1); }
              while ( mpHS->isEmptyPosition(Base::pos()) ); }

    private:
        HashSet const* mpHS;
    };

    explicit HashSet( size_t cap,
                      hash_func_type const& hasher = hash_func_type(),
                      comp_func_type const& comper = comp_func_type(),
                      factory_type const& factory = factory_type() )
    : mHCF(hasher,comper,factory)
    { init(cap); }

    HashSet( HashSet const& that )
    : mHCF(that.mHCF)
    { init(that.size());
      const_iterator end(that.end());
      for ( const_iterator itr(that.begin()); itr != end; ++itr )
        add(itr->begin(),itr->end()); }

    ~HashSet() { destroy(); }

    HashSet& operator=( HashSet const& that )
    { if ( this != &that )
      { clear(); const_iterator end(that.end());
        for ( const_iterator itr(that.begin()); itr != end; ++itr )
          add(itr->begin(),itr->end()); }
      return *this; }

    size_type size() const
    { size_t totSize = 0;
      PPHHS end(mppHHS+mCapacity);
      for ( PPHHS itr(mppHHS); itr != end; ++itr )
        if ( *itr ) totSize += (*itr)->size();
      return totSize; }

    factory_type const& getFactory() const { return mHCF; }

    const_iterator begin() const { return const_iterator(mppHHS,this); }
    const_iterator end() const { return const_iterator(mppHHS+mCapacity,this); }
    const_iterator cbegin() { return const_iterator(mppHHS,this); }
    const_iterator cend() { return const_iterator(mppHHS+mCapacity,this); }
    iterator begin() { return iterator(mppHHS,this); }
    iterator end() { return iterator(mppHHS+mCapacity,this); }

    /// Returns true if the value is in the set.
    /// Not safe to call when there are writers lurking about.  Fine in a
    /// readers-only context.  Sorry about the ridiculous rules for using the
    /// API for this class.
    bool contains( key_type val ) const { return lookup(val); }

    /// Returns pointer to val, or null pointer if val is not in set.
    /// Not safe to call when there are writers lurking about.  Fine in a
    /// readers-only context.  Sorry about the ridiculous rules for using the
    /// API for this class.
    const_pointer lookup( key_type val ) const
    { size_t hash = mHCF.hash(val);
      return findHHSC(hash)->lookup(val,hash); }

    /// Returns true if the value was added (otherwise, it was already present).
    bool add( key_type val )
    { size_t hash = mHCF.hash(val);
      HHS* pHHS = findHHS(hash);
      if ( pHHS->lookup(val)  ) { pHHS->unlock(); return false; }
      insert(val,hash,pHHS); return true; }

    template <class Itr>
    void add( Itr itr, Itr end )
    { while ( itr != end ) add(*itr); }

    /// Adds a novel key quickly.
    /// There is no check for key existence, and there is no protection for
    /// multiple threads.  Use this only when you're sure the key doesn't exist
    /// and you're doing all insertion from a single thread!
    void addDangerously( key_type val )
    { size_t hash = mHCF.hash(val);
      size_t hash2 = hash ^ (hash >> 32);
      size_t mask = mCapacity - 1;
      HHS* pHHS;
      while ( !(pHHS = mppHHS[hash2&mask]) ) mask >>= 1;
      while ( true )
      { try  { pHHS->insertK(val,hash); break; }
        catch ( NoRoomException const& ) { pHHS = split(hash); } } }

    /// Returns a reference to the value.  Inserts value, if not present.
    const_reference operator[]( key_type val )
    { size_t hash = mHCF.hash(val);
      HHS* pHHS = findHHS(hash);
      const_pointer ppp = pHHS->lookup(val,hash);
      if ( ppp ) { pHHS->unlock(); return *ppp; } // EARLY RETURN!
      return insert(val,hash,pHHS); }

    /// Removes the value from the set.  Returns false if value not present.
    bool remove( key_type val )
    { size_t hash = mHCF(val);
      HHS* pHHS = findHHS(hash);
      bool result = pHHS->remove(val,hash);
      pHHS->unlock();
      return result; }

    /// Delete all values.  NB: No locking!  You can't call clear when there
    /// are simultaneous writers.
    HashSet& clear()
    { PPHHS end = mppHHS + mCapacity;
      for ( PPHHS itr(mppHHS); itr != end; ++itr )
        if ( *itr ) (*itr)->clear();
      return *this; }

    bool isEmptyPosition( PPHHS ppHHS ) const
    { size_type idx = ppHHS - mppHHS;
      return idx < mCapacity && !*ppHHS; }

    size_t writeBinary( BinaryWriter& writer ) const
    { size_type cap = mCapacity;
      size_t len = writer.write(cap);
      len += writer.write(mInnerCapacity);
      PPHHS end = mppHHS + mCapacity;
      for ( PPHHS itr = mppHHS; itr != end; ++itr )
      { if ( !*itr ) len += writer.write(false);
        else { len += writer.write(true); len += writer.write(**itr); } }
      return len; }

    void readBinary( BinaryReader& reader )
    { destroy();
      size_type cap; reader.read(&cap);
      ForceAssert(cap && !(cap & (cap-1)));
      mCapacity = cap;
      reader.read(&mInnerCapacity);
      mppHHS = new HHS*[mCapacity];
      memset((void*)mppHHS,0,mCapacity*sizeof(HHS*));
      PPHHS end = mppHHS + mCapacity;
      for ( PPHHS itr = mppHHS; itr != end; ++itr )
      { bool populated; reader.read(&populated);
        if ( populated )
        { *itr = new HHS(mInnerCapacity,mHCF); reader.read(*itr); } } }

    static size_t externalSizeof() { return 0; }

    friend bool operator==( HashSet const& set1, HashSet const& set2 )
    { bool result = set1.size() == set2.size();
      for ( const_iterator oItr(set1.begin()), oEnd(set1.end());
            result && oItr != oEnd; ++oItr )
        for ( ICItr iItr(oItr->begin()), iEnd(oItr->end());
              result && iItr != iEnd; ++iItr )
        { const_pointer pEnt = set2.lookup(*iItr);
          result = pEnt && *pEnt == *iItr; }
      return result; }

    friend bool operator!=( HashSet const& set1, HashSet const& set2 )
    { return !(set1 == set2); }

private:
    void init( size_t cap )
    { int capCeilLg2 = PowerOf2::ceilLg2(4*cap/3);
      int innerCapLg2 = std::max(10,(capCeilLg2+1)/2);
      int outerCapLg2 = 1;
      if ( capCeilLg2 > innerCapLg2 ) outerCapLg2 = capCeilLg2 - innerCapLg2;
      mCapacity = 1ul << outerCapLg2;
      mppHHS = new HHS*[mCapacity];
      PPHHS ppp = mppHHS + mCapacity;
      mInnerCapacity = PowerOf2::getNearbyPrime(innerCapLg2);
      while ( ppp != mppHHS ) *--ppp = new HHS(mInnerCapacity,mHCF);
      mSplitCount = 0; }

    // no locking:  therefore only useful when there are no writers that might
    // cause table reorganizations.
    HHS* findHHSC( size_t hash ) const
    { size_t mask = mCapacity - 1;
      hash ^= hash >> 32;
      HHS* result;
      while ( !(result = mppHHS[hash&mask]) ) mask >>= 1;
      return result; }

    // return with HHS locked
    HHS* findHHS( size_t hash )
    { hash ^= hash >> 32;
      HHS* pHHS;
      while ( true )
      { unsigned splitCount = mSplitCount;
        size_t mask = mCapacity - 1;
        while ( !(pHHS = mppHHS[hash&mask]) ) mask >>= 1;
        pHHS->lock();
        if ( mSplitCount == splitCount ) break;
        pHHS->unlock(); }
      return pHHS; }

    const_reference insert( key_type val, size_t hash, HHS* pHHS )
    { while ( true )
      { try
        { const_reference result = pHHS->insertK(val,hash);
          pHHS->unlock();
          return result; }
        catch ( NoRoomException const& )
        { pHHS = split(hash); } } }

    HHS* split( size_t hash )
    { mHCF.lock();
      size_t mask = mCapacity - 1;
      hash ^= hash >> 32;
      while ( !mppHHS[hash&mask] ) mask >>= 1;
      size_t idxLo = hash&mask++;
      while ( !mppHHS[idxLo|(mask>>1)] ) mask >>= 1;
      size_t idxHi = idxLo | mask;
      if ( idxHi >= mCapacity )
      { PPHHS ppHHS = new HHS*[mCapacity<<1];
        memcpy((void*)ppHHS,(void*)mppHHS,mCapacity*sizeof(HHS*));
        memset((void*)(ppHHS+mCapacity),0,mCapacity*sizeof(HHS*));
        PPHHS ppHHSOld = mppHHS;
        mppHHS = ppHHS; mCapacity <<= 1;
        delete [] ppHHSOld; }
      mHCF.unlock();

      HHS* pOld = mppHHS[idxLo];
      HHS newLo(mInnerCapacity,mHCF);
      HHS* pNewHi = new HHS(mInnerCapacity,mHCF);
      for ( ICItr itr(pOld->begin()), end(pOld->end()); itr != end; ++itr )
      { const_reference val = *itr;
        size_t hash2 = mHCF.hash(val);
        (((hash2^(hash2>>32))&mask)?pNewHi:&newLo)->insertV(val,hash2); }
      swap(*pOld,newLo);
      if ( hash&mask ) pNewHi->lock();
      mppHHS[idxHi] = pNewHi;
      bumpSplitCount();
      if ( hash&mask ) { pOld->unlock(); pOld = pNewHi; }
      return pOld; }

    void bumpSplitCount()
    { unsigned count;
      do
      { count = mSplitCount; }
      while ( !__sync_bool_compare_and_swap(&mSplitCount,count,count+1) ); }

    void destroy()
    { PPHHS end = mppHHS+mCapacity;
      for ( PPHHS itr(mppHHS); itr != end; ++itr )
        delete *itr;
      delete [] mppHHS; }

    PPHHS volatile mppHHS; // length == mCapacity
    size_type volatile mCapacity; // number of HHS*s, not max number of values,
                                  //   always a power of 2
    unsigned mInnerCapacity; // prime number just less than a power of 2
    unsigned volatile mSplitCount; // counter of the number of splits
    HCFLock<T,H,C,F> mHCF; // the spin-lock in this class is used to control
                           // writing to mppHHS
};

template <class T, class H, class C, class F>
struct Serializability< HashSet<T,H,C,F> > : public SelfSerializable {};

#endif /* FEUDAL_HASHSET_H */
