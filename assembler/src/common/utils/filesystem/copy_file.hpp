//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "path_helper.hpp"

namespace fs {

files_t files_by_prefix(std::filesystem::path const& path);
void copy_files_by_prefix(files_t const& files, std::filesystem::path const& to_folder);
void link_files_by_prefix(files_t const& files, std::filesystem::path const& to_folder);
void copy_files_by_ext(std::filesystem::path const& from_folder, std::filesystem::path const& to_folder,
                       std::string const& ext, bool recursive);
}
