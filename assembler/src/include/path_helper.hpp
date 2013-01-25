//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>

#include <string>
#include <vector>

namespace path {

typedef std::vector<std::string> files_t;

bool make_dir(std::string const& folder);
std::string make_temp_dir(const std::string &prefix, const std::string &suffix);
void remove_dir(std::string const& folder);
bool is_regular_file(std::string const& path);
std::string append_path(std::string const& prefix, std::string const& suffix);
std::string current_dir();
void make_full_path(std::string& path);
std::string filename(std::string const& path);
std::string basename(std::string const& path);
std::string extension(std::string const& path);
std::string parent_path(const std::string &path);

// doesn't support symlinks
std::string resolve(std::string const& path);
std::string make_relative_path(const std::string &p, std::string base = current_dir());
}
