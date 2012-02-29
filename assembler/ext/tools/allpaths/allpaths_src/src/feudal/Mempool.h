///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Mempool.h
 * \author tsharpe
 * \date Jul 6, 2009
 *
 * \brief Simple memory sub-allocator.
 */
#ifndef FEUDAL_MEMPOOL_H_
#define FEUDAL_MEMPOOL_H_

#include "system/AlignmentCalculator.h"
#include "system/Assert.h"
#include "system/SpinLockedData.h"
#include <cstddef>
#include <limits>
#include <algorithm>
#include <stack>
#include <vector>

/// Dead-simple memory sub-allocation facility.
/// It hands out chunks of memory, avoiding the overhead associated
/// with tiny allocations, but does so at the cost of not tracking
/// at all any of the memory that's freed.  The Mempool ref-counts
/// its clients, and all the memory in the pool goes away only when
/// all clients have gotten out of the pool.  As stupid as this sounds,
/// it's really exactly what's needed for the majority of our work with
/// double vectors, where the inner vectors, which make oodles of tiny
/// allocations, are loaded and never modified.
class Mempool : public SpinLockedData
{
public:
    Mempool()
    : mpChunk(0), mpPreallocatedChunk(0), mTotalSize(0), mFreeSize(0),
      mChunkSize(DEFAULT_CHUNK_SIZE), mRefCount(0)
    {}

    ~Mempool()
    { if ( mTotalSize ) reportLeak();
      if ( mpPreallocatedChunk )
      { reportUnusedPreallocation(mpPreallocatedChunk);
        killChunkChain(mpPreallocatedChunk); }
      else if ( mpChunk )
        killChunkChain(mpChunk); }

    void* allocate( size_t siz, size_t alignmentReq );
    void free( void* ppp, size_t siz );

    void preAllocate( size_t nBytes, size_t nInstances )
    { if ( !tooBig(nBytes) )
      { nBytes *= nInstances;
        if ( nBytes >= 2*mChunkSize ) preAllocate(nBytes); } }

    size_t ref() { return ++mRefCount; }
    size_t deref() { return --mRefCount; }

    void setChunkSize( size_t chunkSize ) { mChunkSize = chunkSize; }
    size_t getMaxEnchunkableSize() const
    { return mChunkSize/MIN_ALLOCS_PER_CHUNK; }

    bool isUnused() { return mTotalSize == 0 && mRefCount == 0; }

    static size_t const MIN_ALLOCS_PER_CHUNK = 8;

private:
    Mempool( Mempool const& ); // unimplemented -- no copying
    Mempool& operator=( Mempool const& ); // unimplemented -- no copying

    class Chunk
    {
    public:
        Chunk( Chunk* pChunk, size_t siz )
        : mpNext(pChunk),
          mFree(reinterpret_cast<char*>(this+1)),
          mEnd(mFree+siz)
        {}

        // compiler-supplied destructor is OK, since we're killing all chunks
        // via Mempool::killChunkChain

        size_t size() const { return mEnd - start(); }
        size_t freeSize() const { return mEnd - mFree; }

        void* start() { return this+1; }

        void* allocate( size_t siz, size_t alignmentReq )
        { char* result = 0;
          if ( mEnd-siz >= mFree )
          { result = align(alignmentReq);
            mFree = result + siz;
            Assert(mFree <= mEnd); }
          return result; }

        void free( void* addr, size_t siz )
        { char* ppp = reinterpret_cast<char*>(addr);
          if ( ppp+siz == mFree ) mFree = ppp; }

    private:
        Chunk( Chunk const& ); // unimplemented -- no copying
        Chunk& operator=( Chunk const& ); // unimplemented -- no copying

        char const* start() const
        { return reinterpret_cast<char const*>(this+1); }

        char* align( size_t alnReq )
        { Assert(alnReq && !(alnReq & (alnReq-1)) ); // i.e., is a power of 2
          size_t mask = alnReq - 1;
          size_t vvv = (reinterpret_cast<size_t>(mFree)+mask) & ~mask;
          return reinterpret_cast<char*>(vvv); }

        friend class Mempool;
        Chunk* mpNext;
        char* mFree;
        char* mEnd;
    };

    void preAllocate( size_t nBytes );
    bool tooBig( size_t siz ) const { return siz > getMaxEnchunkableSize(); }
    void reportLeak();
    void killChunkChain( Chunk* );
    static void reportUnusedPreallocation( Chunk* );

    Chunk* mpChunk;
    Chunk* mpPreallocatedChunk;
    size_t mTotalSize;
    size_t mFreeSize;
    size_t mChunkSize;
    size_t mRefCount;

    // 2Mb less a little in case exact powers of 2 are inefficient
    static size_t const DEFAULT_CHUNK_SIZE = 2*1024*1024 - 24 - sizeof(Chunk);
};

/// Manages conversion between a short ID and a Mempool.
/// The idea here is that we can't afford the space for a pointer in
/// inner-vector type objects, so we'll just fix the number of pools
/// at 64K and have each outer-vector grab one for use by all its inner-
/// vectors.  We can afford a 16-bit pool ID in the inner-vectors.
/// If you need more than 64K outer-vectors, some of them will share a pool.
class MempoolFinder : SpinLockedData
{
public:
    static MempoolFinder& getInstance() // this class is a singleton
    { return gpInstance ? *gpInstance : createInstance(); }

    unsigned short allocatePool()
    { SpinLocker locker(*this);
      unsigned short result;
      if ( !mFreePools.empty() )
      { result = mFreePools.top(); mFreePools.pop(); }
      else
      { if ( !mNextPool ) mNextPool += 1;
        result = mNextPool++; }
      return result; }

    void freePool( unsigned short poolID )
    { SpinLocker locker(*this);
      mFreePools.push(poolID); }

    Mempool* resolvePool( unsigned short poolID )
    { return mMempools + poolID; }

private:
    MempoolFinder() : mNextPool(1)
    { mMempools[0].setChunkSize(0);
      unsigned short iii = static_cast<unsigned short>(N_POOLS);
      while ( --iii ) mFreePools.push(iii); }
    MempoolFinder( MempoolFinder const& ); // unimplemented -- no copying
    MempoolFinder& operator=( MempoolFinder const& ); // unimplemented -- no copying
    ~MempoolFinder(); // check for leaks

    static MempoolFinder& createInstance() // this class is a singleton
    { static MempoolFinder gInstance;
      gpInstance = &gInstance; return gInstance; }

    static unsigned const N_POOLS = 65536;
    Mempool mMempools[N_POOLS];
    unsigned short mNextPool;
    std::stack<unsigned short,std::vector<unsigned short> > mFreePools;

    static MempoolFinder* gpInstance;
};

template <class T>
class MempoolAllocator
{
public:
    typedef T value_type;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef size_t size_type;
    typedef std::ptrdiff_t difference_type;

    template<typename U>
    struct rebind
    {
        typedef MempoolAllocator<U> other;
    };

    MempoolAllocator() : mPoolID(0) {}
    template<typename U> MempoolAllocator( MempoolAllocator<U> const& that )
    : mPoolID(that.poolID())
    {}

    // compiler-supplied copying and destructor are OK

    pointer address( reference r ) { return &r; }

    const_pointer address( const_reference r ) { return &r; }

    pointer allocate( size_type siz, void* /*hint*/ = 0 )
    { size_t aln = AlignmentCalculator<T>::getAlignment();
      siz *= sizeof(T);
      return reinterpret_cast<pointer>(getPool()->allocate(siz,aln)); }

    void deallocate( pointer p, size_type siz )
    { getPool()->free(p,siz*sizeof(T)); }

    size_type max_size() const
    { return std::numeric_limits<size_type>::max()/sizeof(T); }

    void construct( pointer p, T const& t ) { new (p) T(t); }

    void destroy( pointer p ) { p->~T(); }

    bool operator==( MempoolAllocator const& that )
    { return mPoolID == that.mPoolID; }

    bool operator!=( MempoolAllocator const& that )
    { return mPoolID != that.mPoolID; }

    void preAllocate( size_type nTs, size_type nInstances )
    { if ( nTs && nInstances >= 2*Mempool::MIN_ALLOCS_PER_CHUNK )
        getPool()->preAllocate(nTs*sizeof(T),nInstances); }

    size_t getMaxEnchunkableSize() const
    { return getPool()->getMaxEnchunkableSize()/sizeof(T); }

    unsigned short poolID() const { return mPoolID; }

    friend void swap( MempoolAllocator& alloc1, MempoolAllocator& alloc2 )
    { using std::swap; swap(alloc1.mPoolID,alloc2.mPoolID); }

protected:
    explicit MempoolAllocator( unsigned short poolID ) : mPoolID(poolID) {}

    Mempool* getPool()
    { return finder().resolvePool(mPoolID); }
    Mempool const* getPool() const
    { return finder().resolvePool(mPoolID); }

    static MempoolFinder& finder() { return MempoolFinder::getInstance(); }

private:
    unsigned short mPoolID;
};

template <class T>
class MempoolOwner : public MempoolAllocator<T>
{
    typedef MempoolAllocator<T> Base;
public:
    MempoolOwner()
    : Base(Base::finder().allocatePool())
    { this->getPool()->ref(); }

    MempoolOwner( MempoolOwner const& mo )
    : Base(mo)
    { this->getPool()->ref(); }

    ~MempoolOwner()
    { if ( !this->getPool()->deref() )
        Base::finder().freePool(this->poolID()); }

    friend void swap( MempoolOwner& alloc1, MempoolOwner& alloc2 )
    { swap(static_cast<Base&>(alloc1),static_cast<Base&>(alloc2)); }

private:
    MempoolOwner& operator=( MempoolOwner const& ); // unimplemented -- no assignment
};

#endif /* FEUDAL_MEMPOOL_H_ */
