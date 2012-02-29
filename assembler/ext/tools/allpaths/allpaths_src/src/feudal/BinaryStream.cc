///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file BinaryStream.cc
 * \author tsharpe
 * \date Aug 13, 2009
 *
 * \brief
 */
#include "feudal/BinaryStream.h"
#include "system/System.h"

void BinaryReader::readLoop( char* buf, size_t len )
{
    size_t remain = 0;
    while ( len )
    {
        remain = fillBuf();
        if ( !remain )
            FatalErr("BinaryReader attempted to read past the end of file "
                      << mFR.getFilename());

        if ( remain > len )
            remain = len;
        memcpy(buf, mpBuf, remain);
        buf += remain;
        len -= remain;
    }
    mpBuf += remain;
}

void BinaryReader::testToken()
{
    MagicToken tok;
    if ( !read(&tok).isValid() )
        FatalErr("Reading binary file " << mFR.getFilename()
                  << " failed: Initial token is invalid.");
}
