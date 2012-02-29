///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FileWriter.cc
 * \author tsharpe
 * \date Jan 20, 2012
 *
 * \brief
 */

#include "system/file/FileWriter.h"
#include "system/ErrNo.h"
#include "system/System.h"
#include <sys/select.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

FileWriter const& FileWriter::write( void const* voidbuf, size_t len ) const
{
    char const* buf = static_cast<char const*>(voidbuf);
    size_t nToGo = len;

    while ( nToGo )
    {
        ssize_t nWritten = ::write(mFD,buf,std::min(nToGo,MAX_IO_LEN));
        if ( nWritten == -1 ) // if an error occurred
        {
            ErrNo err;
            if ( err.val() == EAGAIN )
            {
                fd_set fds;
                FD_ZERO(&fds);
                FD_SET(mFD,&fds);
                while ( select(mFD+1,0,&fds,0,0) == -1 )
                {
                    ErrNo err2;
                    if ( err2.val() != EINTR )
                        FatalErr("Select failed on a non-blocking write"
                                    << err2);
                }
            }
            else if ( err.val() == EINTR )
                continue;

            FatalErr("Attempt to write " << len << " bytes to " << mPath <<
                     " failed after writing " << len-nToGo << " bytes" << err);
        }
        nToGo -= nWritten;
        buf += nWritten;
    }

    return *this;
}

int FileWriter::doOpen( char const* path )
{
    int fd = ::open(path,O_RDWR|O_CREAT|O_TRUNC,0664);
    if ( fd == -1 )
    {
        ErrNo err;
        FatalErr("Attempt to open " << path << " for writing failed" << err);
    }
    return fd;
}
