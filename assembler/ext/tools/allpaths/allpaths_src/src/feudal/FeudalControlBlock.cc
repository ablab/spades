///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file FeudalControlBlock.cc
 * \author tsharpe
 * \date Jun 23, 2009
 *
 * \brief
 */
#include "feudal/FeudalControlBlock.h"
#include "system/ErrNo.h"
#include "system/Exit.h"
#include "system/file/FileReader.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <iostream>

using std::cout;
using std::endl;

void FeudalControlBlock::init( FileReader const& fr, bool validate,
                                size_t* pFileLen )
{
    size_t fileLen = fr.getSize();
    if ( pFileLen )
        *pFileLen = fileLen;

    fr.seek(0);
    fr.read(this,sizeof(*this));

    if ( validate && !isValid(fr.getFilename().c_str(),fileLen,true) )
        CRD::exit(1);
}

bool FeudalControlBlock::isValid( char const* fileName, size_t fileLen,
                                    bool verbose ) const
{
    if ( getNFiles() != 1 )
    {
        if ( verbose )
        {
            if ( getNFiles() )
                cout << "Feudal file " << fileName << " is in " << getNFiles()
                     << " files, but we require it to be in a single file."
                     << endl;
            else
                cout << "Feudal file " << fileName << " claims to exist in 0 "
                        "files.  It's probably not a feudal file, or the "
                        "process of creating it was interrupted." << endl;
        }
        return false;
    }

    size_t offsetBytes = mFixedOffset - mVarOffset;
    if ( offsetBytes % sizeof(size_t) )
    {
        if ( verbose )
            cout << "Feudal file " << fileName
                 << " has offset info that does not contain an integral number "
                    "of offsets." << endl;
        return false;
    }

    //unsigned int nnn = offsetBytes/sizeof(size_t) - 1UL;
    if ( ((offsetBytes/sizeof(size_t) - 1UL) & 0xffffffffUL) != mN )
    {
        if ( verbose )
            cout << "Feudal file " << fileName
                 << " doesn't have the right number of offsets for the number "
                    "of elements." << endl;
        return false;
    }

    if ( !isCompressed() )
    {
        size_t fixedBytes = fileLen - mFixedOffset;
        size_t nElements = getNElements();
        if ( !nElements )
        {
            if ( fixedBytes )
            {
                if ( verbose )
                    cout << "Empty feudal file " << fileName
                         << " has some fixed data.  Must be garbage at the "
                            "file's end." << endl;
                return false;
            }
        }
        else
        {
            if ( fixedBytes % nElements )
            {
                if ( verbose )
                    cout << "Feudal file " << fileName
                         << " has fixed info that is not an integral number of "
                            "bytes per element." << endl;
                return false;
            }

            if ( mSizeofFixed && ((fixedBytes/nElements)&0xff) != mSizeofFixed )
            {
                if ( verbose )
                    cout << "Feudal file " << fileName
                         << " doesn't have the right amount of fixed data for "
                            "the number of elements." << endl;
                return false;
            }
        }
    }

    return true;
}

bool FeudalControlBlock::isGoodFeudalFile( char const* fileName, bool verbose )
{
    size_t fileLen;
    FeudalControlBlock fcb(fileName,false,&fileLen);
    return fcb.isValid(fileName,fileLen,verbose);
}
