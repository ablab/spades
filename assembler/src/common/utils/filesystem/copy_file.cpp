//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "copy_file.hpp"

#include "utils/logger/logger.hpp"

#include <boost/algorithm/string.hpp>

#include <iostream>
#include <fstream>

#include <unistd.h>
#include <dirent.h>

#include <sys/stat.h>
#include <sys/types.h>

namespace fs {

namespace details {

void hard_link(std::filesystem::path from_path, std::filesystem::path to_path) {
    from_path = make_full_path(from_path);
    to_path = make_full_path(to_path);

    if (from_path == to_path)
        return;

    if (link(from_path.c_str(), to_path.c_str()) == -1) {
        WARN("Failed to create link. Reason: " << strerror(errno) << ". Error code: " << errno << ". Copying instead");
        copy_file(from_path, to_path);
    }
}

files_t files_in_folder(std::filesystem::path const& path) {
    DIR *dp;
    if ((dp = opendir(path.c_str())) == NULL)
        throw std::runtime_error("can not open folder " + path.native());

    files_t files;

    struct dirent *dirp;
    while ((dirp = readdir(dp)) != NULL)
        if (dirp->d_type == DT_REG)
            files.push_back(path / dirp->d_name);

    closedir(dp);
    return files;
}

files_t folders_in_folder(std::filesystem::path const& path) {
    DIR *dp;
    if ((dp  = opendir(path.c_str())) == NULL)
        throw std::runtime_error("can not open folder " + path.native());

    files_t folders;

    struct dirent *dirp;
    while ((dirp = readdir(dp)) != NULL)
        if (dirp->d_type == DT_DIR) {
            std::filesystem::path folder = dirp->d_name;

            if (folder != "." && folder != "..")
                folders.push_back(path / folder);
        }

    closedir(dp);
    return folders;
}

} // details

files_t files_by_prefix(std::filesystem::path const& path) {
    using namespace details;
    files_t files;

    std::filesystem::path folder(path.parent_path());
    std::filesystem::path prefix = path.filename();

    files_t out_files;
    const files_t all_files = files_in_folder(folder);

    for (auto it = all_files.begin(); it != all_files.end(); ++it) // no std::copy_if before C++11
        if (boost::starts_with(it->filename().c_str(), prefix.c_str()))
            out_files.push_back(*it);

    return out_files;
}

void copy_files_by_prefix(files_t const& files, std::filesystem::path const& to_folder) {
    using namespace details;

    for (auto it = files.begin(); it != files.end(); ++it) {
        files_t files_to_copy = files_by_prefix(*it);

        for (auto it = files_to_copy.begin(); it != files_to_copy.end(); ++it)
            copy_file(*it, to_folder / it->filename());
    }
}

void link_files_by_prefix(files_t const& files, std::filesystem::path const& to_folder) {
    using namespace details;

    for (auto it = files.begin(); it != files.end(); ++it) {
        files_t files_to_copy = files_by_prefix(*it);

        for (auto it = files_to_copy.begin(); it != files_to_copy.end(); ++it)
            hard_link(*it, to_folder / it->filename());
    }
}

void copy_files_by_ext(std::filesystem::path const& from_folder,
                       std::filesystem::path const& to_folder,
                       std::string const& ext, bool recursive) {
    using namespace details;

    files_t files = files_in_folder(from_folder);

    for (auto it = files.begin(); it != files.end(); ++it)
        if (it->extension() == ext)
            copy_file(*it, to_folder / it->filename());

    if (recursive) {
        files_t folders = folders_in_folder(from_folder);

        for (auto it = folders.begin(); it != folders.end(); ++it) {
            std::filesystem::path subdir = to_folder / it->filename();
            create_directory(subdir);
            copy_files_by_ext(*it, subdir, ext, recursive);
        }
    }
}

}
