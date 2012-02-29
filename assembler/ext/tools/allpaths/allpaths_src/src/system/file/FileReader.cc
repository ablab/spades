///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FileReader.cc
 * \author tsharpe
 * \date Jan 12, 2012
 *
 * \brief
 */

#include "system/file/FileReader.h"
#include "system/ErrNo.h"
#include "system/System.h"
#include <algorithm>
#include <cerrno>
#include <sys/mman.h>
#include <sys/select.h>
#include <sys/types.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

size_t const FileReader::MAX_IO_LEN;

size_t FileReader::readOnce( void* voidbuf, size_t len ) const
{
    char* buf = static_cast<char*>(voidbuf);
    ssize_t nRead;
    while ( (nRead = ::read(mFD,buf,std::min(len,MAX_IO_LEN))) == -1 )
    {
        ErrNo err;
        if ( err.val() == EAGAIN )
            wait();
        else if ( err.val() != EINTR )
            FatalErr("Attempt to read " << len << " bytes from " << mPath <<
                        " failed" << err);
    }
    return nRead;
}

size_t FileReader::readSome( void* voidbuf, size_t len ) const
{
    char* buf = static_cast<char*>(voidbuf);
    size_t nToGo = len;

    while ( nToGo )
    {
        ssize_t nRead = ::read(mFD,buf,std::min(nToGo,MAX_IO_LEN));
        if ( !nRead ) // at EOF
            break;

        if ( nRead == -1 ) // if an error occurred
        {
            ErrNo err;
            if ( err.val() == EAGAIN )
                wait();
            else if ( err.val() == EINTR )
                continue;

            FatalErr("Attempt to read " << len << " bytes from " << mPath <<
                     " failed after reading " << len-nToGo << " bytes" << err);
        }
        nToGo -= nRead;
        buf += nRead;
    }

    return len-nToGo;
}

FileReader const& FileReader::read( void* buf, size_t len ) const
{
    size_t nRead = readSome(buf,len);
    if ( nRead < len )
        FatalErr("Attempt to read " << len << " bytes from " << mPath
                << " failed.  There were " << nRead << " bytes before EOF.");
    return *this;
}

FileReader const& FileReader::seek( size_t off ) const
{
    if ( ::lseek(mFD,off,SEEK_SET) == -1 )
    {
        ErrNo err;
        FatalErr("Attempt to lseek " << mPath << " failed" << err);
    }
    return *this;
}

FileReader const& FileReader::seekRel( long off ) const
{
    if ( ::lseek(mFD,off,SEEK_CUR) == -1 )
    {
        ErrNo err;
        FatalErr("Attempt to lseek(SEEK_CUR) " << mPath << " failed" << err);
    }
    return *this;
}

FileReader const& FileReader::seekEnd( long off ) const
{
    if ( ::lseek(mFD,off,SEEK_END) == -1 )
    {
        ErrNo err;
        FatalErr("Attempt to lseek(SEEK_END) " << mPath << " failed" << err);
    }
    return *this;
}

struct stat FileReader::getStat() const
{
    struct stat sb;
    if ( fstat(mFD,&sb) == -1 )
    {
        ErrNo err;
        FatalErr("Can't fstat " << mPath << err );
    }
    return sb;
}

void* FileReader::map( size_t offset, size_t len, bool readOnly )
{
    int prot = PROT_READ;
    if ( !readOnly )
        prot |= PROT_WRITE;
    void* addr = mmap(0,len,prot,MAP_SHARED,mFD,offset);
    if ( addr == MAP_FAILED )
    {
        ErrNo err;
        FatalErr("Unable to mmap " << mPath << err);
    }
    return addr;
}

void FileReader::doOpen()
{
    if ( (mFD = ::open(mPath.c_str(),O_RDONLY)) == -1 )
    {
        ErrNo err;
        FatalErr("Attempt to open " << mPath << " for reading failed" << err);
    }
    mMyFD = true;
}

void FileReader::doClose()
{
    while ( ::close(mFD) == -1 )
    {
        ErrNo err;
        if ( err.val() != EINTR )
            FatalErr("Attempt to close " << mPath << " failed" << err);
    }
    mFD = -1;
}

void FileReader::wait() const
{
    fd_set fds;
    FD_ZERO(&fds);
    FD_SET(mFD,&fds);
    while ( select(mFD+1,&fds,0,0,0) == -1 )
    {
        ErrNo err;
        if ( err.val() != EINTR )
            FatalErr("Select failed on a non-blocking read" << err);
    }
}
