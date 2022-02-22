//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <string>
#include <vector>
#include <filesystem>

namespace fs {
//todo review and make names consistent!

typedef std::vector<std::filesystem::path> files_t;

std::filesystem::path make_temp_dir(std::filesystem::path const &prefix, std::string const &suffix);

std::filesystem::path make_full_path(std::filesystem::path const &path);

std::filesystem::path screen_whitespaces(std::filesystem::path const &path);

std::filesystem::path resolve(const std::filesystem::path &path);

std::filesystem::path make_launch_time_dir_name();
}
