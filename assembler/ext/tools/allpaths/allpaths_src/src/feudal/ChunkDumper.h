///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file ChunkDumper.h
 * \author tsharpe
 * \date Sep 22, 2010
 *
 * \brief
 */
#ifndef FEUDAL_CHUNKDUMPER_H_
#define FEUDAL_CHUNKDUMPER_H_

#include "feudal/BinaryStream.h"
#include "system/Assert.h"
#include "system/SpinLockedData.h"
#include <cstddef>
#include <cstring>
#include <iostream>

/// Write binary file in chunks.
/// This stupid thing just buffers chunks of stuff that we're writing
/// out to a binary file, making sure that the chunks get written in the
/// correct order.  It's (vaguely) useful in multi-threaded contexts when
/// you're producing a vector of a known size, but each thread is producing
/// the data a chunk at a time.  In this case, you can save memory by writing
/// the data chunk-wise, rather than accumulating the whole mess and then
/// writing it.
/// The type V must be vector-like.  It must be default constructable,
/// swappable, have size_type and value_type typedefs, have a size() method,
/// a reserve() method, and have an operator[] that returns a legitimate
/// reference to an element.
template <class V>
class ChunkDumper : public SpinLockedData
{
public:
    typedef typename V::size_type size_type;
    typedef typename V::value_type value_type;

    ChunkDumper( char const* filename, size_type nValsTotal, size_t nChunks )
    : mWriter(filename), mNValsTotal(nValsTotal), mNChunks(nChunks),
      mChunks(new V*[nChunks]), mNextChunk(0), mNValsWritten(0)
    { memset(mChunks,0,nChunks*sizeof(V*)); mWriter.write(nValsTotal); }

    ~ChunkDumper()
    { if ( mChunks ) close(); }

    void close()
    { ForceAssertEq(mNChunks,mNextChunk);
      ForceAssertEq(mNValsWritten,mNValsTotal);
      mWriter.close();
      delete [] mChunks;
      mChunks = 0; }

    /// Stores your vals for writing, and swaps in an empty V in its place.
    void dumpChunk( size_t chunk, V& vals )
    { ForceAssertGe(chunk,mNextChunk);
      V* pV = new V();
      using std::swap; swap(*pV,vals);
      SpinLocker lock(*this);
      mChunks[chunk] = pV;
      if ( chunk == mNextChunk ) doDump(); }

private:
    ChunkDumper( ChunkDumper const& ); // unimplemented -- no copying
    ChunkDumper& operator=( ChunkDumper const& ); // unimplemented -- no copying

    void doDump()
    { while ( mNextChunk < mNChunks )
      { V* pV = mChunks[mNextChunk];
        if ( !pV ) break;
        if ( pV->size() )
        { value_type const* start = &(*pV)[0];
          value_type const* end = start + pV->size();
          mWriter.write(start,end);
          mNValsWritten += pV->size(); }
        delete pV; mChunks[mNextChunk++] = 0; } }

    BinaryWriter mWriter;
    size_t mNValsTotal;
    size_t mNChunks;
    V** mChunks;
    size_t mNextChunk;
    size_t mNValsWritten;
};

#endif /* FEUDAL_CHUNKDUMPER_H_ */
