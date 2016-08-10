//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/path_helper.hpp"
#include <string>

namespace path {

path::files_t files_by_prefix(std::string const& path);
void copy_files_by_prefix(path::files_t const& files, std::string const& to_folder);
void link_files_by_prefix(path::files_t const& files, std::string const& to_folder);
void copy_files_by_ext(std::string const& from_folder, std::string const& to_folder, std::string const& ext, bool recursive);

}
