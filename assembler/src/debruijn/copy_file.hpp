//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * copy_file.hpp
 *
 *  Created on: 8 Sep 2011
 *      Author: valery
 */

#pragma once

#include <unistd.h>
#include <dirent.h>

#include <sys/stat.h>
#include <sys/types.h>

#include "path_helper.hpp"

namespace details
{

using namespace path;

inline void copy_file(string from_path, string to_path)
{
    using namespace std;

    make_full_path(from_path);
    make_full_path(to_path  );

    ifstream source(from_path, ios::binary);
    ofstream dest  (to_path  , ios::binary);

    dest << source.rdbuf();
}


inline void hard_link(string from_path, string to_path)
{
    make_full_path(from_path);
    make_full_path(to_path  );

    if (link(from_path.c_str(), to_path.c_str()) != 0)
        copy_file(from_path, to_path);
}

inline files_t files_in_folder(string const& path)
{
    DIR *dp;
    if ((dp = opendir(path.c_str())) == NULL)
        throw std::runtime_error("can not open folder " + path);

    files_t files;

    struct dirent *dirp;
    while ((dirp = readdir(dp)) != NULL)
        if (dirp->d_type == DT_REG)
            files.push_back(append_path(path, dirp->d_name));

    closedir(dp);
    return files;
}

inline files_t folders_in_folder(string const& path)
{
    DIR *dp;
    if ((dp  = opendir(path.c_str())) == NULL)
        throw std::runtime_error("can not open folder " + path);

    files_t folders;

    struct dirent *dirp;
    while ((dirp = readdir(dp)) != NULL)
        if (dirp->d_type == DT_DIR)
        {
            string folder = dirp->d_name;

            if (folder != "." &&
                folder != "..")
            {
                folders.push_back(append_path(path, folder));
            }
        }

    closedir(dp);
    return folders;
}

} // details

inline details::files_t files_by_prefix(string const& path)
{
    using namespace details;
    files_t files;

    string folder(parent_path(path));
    std::string prefix = filename(path);

          files_t out_files;
    const files_t all_files = files_in_folder(folder);

    for (auto it = all_files.begin(); it != all_files.end(); ++it) // no std::copy_if before C++11
        if (boost::starts_with(filename(*it), prefix))
            out_files.push_back(*it);

    return out_files;
}

inline void copy_files_by_prefix(details::files_t const& files, string const& to_folder)
{
    using namespace details;

    for (auto it = files.begin(); it != files.end(); ++it)
    {
        files_t files_to_copy = files_by_prefix(*it);

        for (auto it = files_to_copy.begin(); it != files_to_copy.end(); ++it)
            copy_file(*it, append_path(to_folder, filename(*it)));
    }
}

inline void link_files_by_prefix(details::files_t const& files, string const& to_folder)
{
    using namespace details;

    for (auto it = files.begin(); it != files.end(); ++it)
    {
        files_t files_to_copy = files_by_prefix(*it);

        for (auto it = files_to_copy.begin(); it != files_to_copy.end(); ++it)
            hard_link(*it, append_path(to_folder, filename(*it)));
    }
}

inline void copy_files_by_ext(string const& from_folder, string const& to_folder, std::string const& ext, bool recursive)
{
    using namespace details;

    files_t files = files_in_folder(from_folder);

    for (auto it = files.begin(); it != files.end(); ++it)
        if (boost::ends_with(*it, ext))
            copy_file(*it, append_path(to_folder, filename(*it)));

    if (recursive)
    {
        files_t folders = folders_in_folder(from_folder);

        for (auto it = folders.begin(); it != folders.end(); ++it)
        {
            string subdir = append_path(to_folder, filename(*it));
            path:: make_dir(subdir);
            copy_files_by_ext(*it, subdir, ext, recursive);
        }
    }
}
