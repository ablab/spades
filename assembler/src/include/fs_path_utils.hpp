//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * fs_path_utils.hpp
 *
 *  Created on: Jan 25, 2012
 *      Author: valery
 */

#pragma once

namespace boost
{
namespace filesystem
{

// copied from http://stackoverflow.com/questions/1746136/how-do-i-normalize-a-pathname-using-boostfilesystem
inline fs::path resolve(const fs::path& p)
{
    using namespace fs;

    path result;
    for(path::iterator it= p.begin(); it != p.end(); ++it)
    {
        if(*it == "..")
        {
            // /a/b/.. is not necessarily /a if b is a symbolic link
            if(boost::filesystem::is_symlink(result) )
                result /= *it;
            // /a/b/../.. is not /a/b/.. under most circumstances
            // We can end up with ..s in our result because of symbolic links
            else if(result.filename() == "..")
                result /= *it;
            // Otherwise it should be safe to resolve the parent
            else
                result = result.parent_path();
        }
        else if(*it == ".")
        {
            // Ignore
        }
        else
        {
            // Just cat other path entries
            result /= *it;
        }
    }
    return result;
}

inline fs::path make_relative_path(fs::path p, fs::path base = fs::initial_path())
{
    p    = resolve(p);
    base = resolve(base);

    fs::path pp = p.parent_path();

    auto i = pp  .begin();
    auto j = base.begin();

    while (i != pp.end() && j != base.end() && *i == *j)
    {
        ++i;
        ++j;
    }

    fs::path result;
    for ( ; j != base.end(); ++j)
        result = fs::path("..") / result;

    for ( ; i != pp.end(); ++i)
        result /= *i;

    return result / p.filename();
}

} // filesystem
} // boost
