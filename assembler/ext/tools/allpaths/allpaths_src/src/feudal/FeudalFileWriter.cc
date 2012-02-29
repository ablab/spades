///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file FeudalFileWriter.cc
 * \author tsharpe
 * \date Mar 13, 2009
 *
 * \brief Utility for writing Feudal Files incrementally.
 */
#include "feudal/FeudalFileWriter.h"
#include "system/Exit.h"
#include <iostream>
using std::cout;
using std::endl;

FeudalControlBlock FeudalFileWriter::gInvalidFCB(0,0,0,0,0,0);

FeudalFileWriter::FeudalFileWriter( char const* filename,
                                    FeudalFileWriter::size_type vecSize,
                                    FeudalFileWriter::size_type eltSize,
                                    FeudalFileWriter::size_type fixedLenDataLen,
                                    unsigned long estimatedNElements )
: mWriter(filename,false),
  mNextOffset(sizeof(FeudalControlBlock)),
  mVecSize(vecSize),
  mEltSize(eltSize),
  mFixedLenDataLen(fixedLenDataLen)
{
    mWriter.write(gInvalidFCB);
    mOffsets.reserve(estimatedNElements);
    mFixedLenData.reserve(fixedLenDataLen*estimatedNElements);
}

FeudalFileWriter::~FeudalFileWriter()
{
    if ( mWriter.isOpen() )
        close();
}

void FeudalFileWriter::addElement( size_t varLen, void const* fixedLenData )
{
    mOffsets.push_back(mNextOffset);
    mNextOffset += varLen;

    if ( mFixedLenDataLen )
    {
        char const* fStart = reinterpret_cast<char const*>(fixedLenData);
        mFixedLenData.insert( mFixedLenData.end(), fStart,
                                fStart+mFixedLenDataLen );
    }
}

void FeudalFileWriter::checkPoint()
{
    if ( !mWriter.isOpen() )
    {
        cout << "You can't checkpoint a FeudalFileWriter after you've closed "
                "it." << endl;
        CRD::exit(1);
    }
    finish();
    mWriter.seek(mNextOffset);
}

void FeudalFileWriter::close()
{
    if ( !mWriter.isOpen() )
        cout << "Warning:  closing feudal file " << mWriter.getFilename()
             << " that is already closed." << endl;
    else
    {
        finish();
        mWriter.close();
    }
}

void FeudalFileWriter::finish()
{
    if ( mOffsets.size() )
        mWriter.write(&*mOffsets.begin(),&*mOffsets.end());
    mWriter.write(mNextOffset);
    if ( mFixedLenData.size() )
        mWriter.write(&*mFixedLenData.begin(),&*mFixedLenData.end());

    mWriter.flush();
    mWriter.seek(0ul);
    FeudalControlBlock fcb(mOffsets.size(),
                           mNextOffset-sizeof(FeudalControlBlock),
                           mFixedLenDataLen, mVecSize, mEltSize);
    mWriter.write(fcb);
}
