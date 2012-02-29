///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file Directory.cc
 * \author tsharpe
 * \date Dec 5, 2011
 */
#include "system/file/Directory.h"
#include "system/ErrNo.h"
#include "system/System.h"
#include <cerrno>
#include <sys/stat.h>
#include <unistd.h>

#define ITR Directory::const_iterator

ITR::const_iterator( std::string const& path )
: mDir(path), mpDIR(opendir(mDir.c_str())), mPos(-1L)
{
    checkStream();
    nextEntry();
}

ITR& ITR::operator=( const_iterator const& that )
{
    endStream();
    mDir = that.mDir;
    if ( that.mpDIR )
    {
        mpDIR = opendir(mDir.c_str());
        checkStream();
        seekdir(mpDIR,that.mPos);
        nextEntry();
    }
    return *this;
}

void ITR::nextEntry()
{
    mPos = telldir(mpDIR);
    if ( mPos == -1 )
    {
        ErrNo err;
        FatalErr("Can't get directory position for " << mDir << err);
    }

    errno = 0;
    struct dirent* pDirEnt = readdir(mpDIR);
    if ( pDirEnt ) // got an entry
        mFile = File(mDir+pDirEnt->d_name);
    else if ( !errno ) // at EOF
        endStream();
    else // error
    {
        ErrNo err;
        FatalErr("Can't read directory " << mDir << err);
    }
}

void ITR::checkStream()
{
    if ( !mpDIR )
    {
        ErrNo err;
        FatalErr("Can't open directory " << mDir << err);
    }
}

void ITR::endStream()
{
    if ( mpDIR )
    {
        if ( closedir(mpDIR) )
        {
            ErrNo err;
            FatalErr("Can't close directory " << mDir << err);
        }
        mpDIR = 0;
    }
    mPos = -1L;
}

bool Directory::create( bool recursive, int mode ) const
{
    if ( isValid() )
        return false;

    if ( recursive )
        directory().create(true,mode);

    if ( mkdir(toString().c_str(),mode) )
    {
        ErrNo err;
        FatalErr("Can't create directory " << toString() << err);
    }
    clearStat();
    return true;
}

void Directory::remove() const
{
    if ( rmdir(toString().c_str()) == -1 )
    {
        ErrNo err;
        FatalErr("Can't remove file " << mPath << err);
    }
    clearStat();
}
