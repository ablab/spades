/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/**
 * \file File.cc
 * \author aaron
 * \date Wednesday, December 03, 2008
 *
 * \brief our basic file class, that describes a system file
 *
 * our basic file class, that describes a system file
 *
 */

#include "system/file/File.h"
#include "system/file/Directory.h"
#include "system/file/SymLink.h"
#include "system/ErrNo.h"
#include "system/System.h"
#include <unistd.h>

std::string File::CWD(".");
std::string File::ROOT("/");

Directory File::asDir() const
{
    return Directory(mPath);
}

SymLink File::asLink() const
{
    return SymLink(mPath);
}

Directory File::directory() const
{
    return Directory(dirname());
}

void File::remove() const
{
    if ( ::remove(mPath.c_str()) == -1 )
    {
        ErrNo err;
        FatalErr("Can't remove file " << mPath << err);
    }
    clearStat();
}

void File::rename( File const& file ) const
{
    if ( ::rename(mPath.c_str(),file.mPath.c_str()) == -1 )
    {
        ErrNo err;
        FatalErr("Can't rename file " << mPath << " to " << file.mPath << err);
    }
    clearStat();
    file.clearStat();
}

void File::patchPath()
{
    size_t pos = 0;
    while ( (pos = mPath.find("//",pos)) != std::string::npos )
        mPath.erase(pos,1);
    pos = 0;
    while ( (pos = mPath.find("/./",pos)) != std::string::npos )
        mPath.erase(pos,2);
    pos = mPath.size();
    if ( pos > 1 && mPath.find("/.",pos-2) == pos-2 )
        mPath.resize(pos-1);
    pos = mPath.size();
    if ( pos > 1 && mPath[pos-1]=='/' )
        mPath.resize(pos-1);
}

void File::setStat() const
{
    if ( !mpStat ) mpStat = new struct stat;
    if ( lstat(mPath.c_str(),mpStat) == -1 )
    {
        mType = NOT_A_FILE;
        delete mpStat;
        mpStat = 0;
    }
    else
    {
        mType = static_cast<FILETYPE>(mpStat->st_mode & S_IFMT);
        if ( mType == SYM_LINK && ::stat(mPath.c_str(),mpStat) == -1 )
        {
            delete mpStat;
            mpStat = 0;
        }
    }
}
