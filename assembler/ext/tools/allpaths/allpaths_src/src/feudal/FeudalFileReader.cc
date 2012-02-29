///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
 * \file FeudalFileReader.cc
 * \author tsharpe
 * \date Wednesday, December 03, 2008
 *
 * \brief FeudalFileReader implementation
 */
#include "feudal/FeudalFileReader.h"
#include "system/Assert.h"
#include "system/ErrNo.h"
#include "system/file/FileReader.h"
#include "system/System.h"
#include <unistd.h>
#include <sys/mman.h>

FeudalFileReader::FeudalFileReader( char const* filename )
: mpMapper(new Mapper(filename)), mCurEle(mpMapper->getNElements()),
  mReader(filename,false)
{
    mpMapper->ref();
}

FeudalFileReader::FeudalFileReader( FeudalFileReader const& that )
: mpMapper(that.mpMapper), mCurEle(mpMapper->getNElements()),
  mReader(that.mReader.getFilename().c_str(),false)
{
    SpinLocker locker(*mpMapper);
    mpMapper->ref();
}

FeudalFileReader::~FeudalFileReader()
{
    size_t count;
    if ( true )
    {
        SpinLocker locker(*mpMapper);
        count = mpMapper->deref();
    }
    if ( !count )
        delete mpMapper;
}

size_t FeudalFileReader::getPreallocation( size_t maxSize, size_t extSize,
                                           size_t start, size_t end ) const
{
    AssertLe(end,getNElements());
    AssertLe(start,end);

    size_t result = 0;
    while ( start < end )
    {
        size_t nnn = getDataLen(start++)/extSize;
        if ( nnn <= maxSize )
            result += nnn;
    }
    return result;
}

FeudalFileReader::Mapper::Mapper( char const* filename )
: mFilename(filename), mRefCount(0), mFCB(0,0,0,0,0)
{
    ForceAssertEq(sizeof(size_t),8u);

    FileReader fr(filename);

    size_t fileLen;
    new (&mFCB) FeudalControlBlock(fr,true,&fileLen);

    size_t varTabOffset = mFCB.getVarTabOffset();
    size_t mappedOffset = varTabOffset - varTabOffset%getpagesize();
    mMappedLen = fileLen - mappedOffset;
    mMappedBit = static_cast<char*>(fr.map(mappedOffset,mMappedLen,true));
    mOffsetsTable = mMappedBit + (varTabOffset - mappedOffset);
    mFixedData = mMappedBit + (mFCB.getFixedOffset() - mappedOffset);
}

FeudalFileReader::Mapper::~Mapper()
{
    munmap(mMappedBit,mMappedLen);
}
