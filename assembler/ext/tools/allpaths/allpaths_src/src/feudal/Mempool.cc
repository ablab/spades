///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Mempool.cc
 * \author tsharpe
 * \date Jul 6, 2009
 *
 * \brief Simple memory sub-allocator.
 */
#include "feudal/Mempool.h"
#include <iostream>

using std::cout;
using std::endl;

MempoolFinder* MempoolFinder::gpInstance;

void* Mempool::allocate( size_t siz, size_t alignmentReq )
{
    if ( tooBig(siz) )
        return new char[siz];

    SpinLocker locker(*this);
    void* result;
    if ( !mpChunk || !(result = mpChunk->allocate(siz,alignmentReq)) )
    {
        if ( mpPreallocatedChunk )
        {
            mpChunk = mpPreallocatedChunk;
            mpPreallocatedChunk = 0;
        }
        else
        {
            void* ppp = new char[mChunkSize+sizeof(Chunk)];
            mpChunk = new (ppp) Chunk(mpChunk,mChunkSize);
            mTotalSize += mChunkSize;
            mFreeSize += mChunkSize;
        }
        result = mpChunk->allocate(siz,alignmentReq);
    }

    mFreeSize -= siz;
    AssertLe(mFreeSize,mTotalSize);
    return result;
}

void Mempool::free( void* ppp, size_t siz )
{
    if ( tooBig(siz) )
    {
        delete [] static_cast<char*>(ppp);
        return;
    }

    Chunk* pPre = 0;
    Chunk* pChunk = 0;

    if ( true )
    {
        SpinLocker locker(*this);
        mFreeSize += siz;
        mpChunk->free(ppp,siz);
        AssertLe(mFreeSize,mTotalSize);
        if ( mFreeSize >= mTotalSize )
        {
            pPre = mpPreallocatedChunk;
            pChunk = mpChunk;
            mpPreallocatedChunk = mpChunk = 0;
            mTotalSize = mFreeSize = 0;
        }
    }

    if ( pPre )
    {
        reportUnusedPreallocation(pPre);
        killChunkChain(pPre);
    }
    else if ( pChunk )
    {
        killChunkChain(pChunk);
    }
}

void Mempool::preAllocate( size_t nBytes )
{
    char* ppp = new char[nBytes+sizeof(Chunk)];

    SpinLocker locker(*this);
    if ( mpPreallocatedChunk )
        delete [] ppp;
    else
    {
        mpPreallocatedChunk = new (ppp) Chunk(mpChunk,nBytes);
        mTotalSize += nBytes;
        mFreeSize += nBytes;
    }
}

void Mempool::reportLeak()
{
    cout << "Warning: Mempool is reclaiming " <<
            mTotalSize-mFreeSize << " leaked bytes." << endl;
}

void Mempool::killChunkChain( Chunk* pChunk )
{
    if ( pChunk->mpNext )
        killChunkChain(pChunk->mpNext);
    pChunk->~Chunk();
    delete [] reinterpret_cast<char*>(pChunk);
}

void Mempool::reportUnusedPreallocation( Chunk* pChunk )
{
    cout << "Warning: Mempool has an unused " << pChunk->size() <<
            "-byte preallocation." << endl;
}

MempoolFinder::~MempoolFinder()
{
    unsigned int trouble = 0;
    Mempool* end = mMempools + N_POOLS;
    for ( Mempool* pPool = mMempools; pPool != end; ++pPool )
    {
        if ( pPool->ref() > 1 )
            trouble += 1;
    }
    if ( trouble )
        cout << "Warning: There were " << trouble <<
        " leaked mempools which were cleaned up at the end of the program run."
        << endl;
}
