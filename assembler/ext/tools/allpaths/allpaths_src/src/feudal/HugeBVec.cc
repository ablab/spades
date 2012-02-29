///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file HugeBVec.cc
 * \author tsharpe
 * \date Sep 16, 2010
 *
 * \brief
 */
#include "feudal/HugeBVec.h"
#include "system/ErrNo.h"
#include "system/file/FileReader.h"
#include "system/System.h"
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

HugeBVec::HugeBVec( char const* fileName )
{
    FileReader fr(fileName);
    size_t fileLen = fr.getSize();
    void const* addr = fr.map(0,fileLen,true);
    mpBuf = static_cast<unsigned char const*>(addr);
    mSize = *static_cast<size_type const*>(addr);
    mpBuf += sizeof(mSize);

    if ( fileLen != physSize()+sizeof(mSize) )
        FatalErr("File size doesn't match recorded size in HugeBVec " <<
                fileName);
}

HugeBVec::~HugeBVec()
{
    void* addr = const_cast<unsigned char*>(mpBuf-sizeof(mSize));
    munmap(addr,physSize()+sizeof(mSize));
}
