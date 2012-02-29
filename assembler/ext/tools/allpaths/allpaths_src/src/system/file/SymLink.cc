/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file SymLink.cc
 * \author tsharpe
 * \date Feb 17, 2009
 *
 * \brief
 *
 *
 */
#include "system/file/SymLink.h"
#include "system/ErrNo.h"
#include "system/System.h"
#include <unistd.h>

File SymLink::target() const
{
    size_t const bufSiz = 8193;
    char buf[bufSiz];
    int len = readlink(toString().c_str(),buf,bufSiz-1);
    if ( len < 0 )
    {
        ErrNo err;
        FatalErr("Can't read symlink " << toString() << err);
    }
    buf[len] = 0;
    return File(buf);
}

void SymLink::setTarget( File const& target, bool force ) const
{
    char const* fileName = toString().c_str();
    if ( force && isLink() && unlink(fileName) == -1 )
    {
        ErrNo err;
        FatalErr("Can't remove existing symlink " << toString()
                    << " in order to retarget it" << err);
    }
    if ( symlink(target.toString().c_str(), fileName) == -1 )
    {
        ErrNo err;
        FatalErr("Can't symlink " << toString() << err);
    }
    clearStat();
}
