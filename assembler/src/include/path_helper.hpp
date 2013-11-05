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
#include "verify.hpp"

namespace path {

typedef std::vector<std::string> files_t;

bool make_dir(std::string const& folder);
std::string make_temp_dir(std::string const& prefix, std::string const& suffix);
void remove_dir(std::string const& folder);
bool is_regular_file(std::string const& path);
std::string append_path(std::string const& prefix, std::string const& suffix);
std::string current_dir();
void make_full_path(std::string& path);
std::string filename(std::string const& path);
std::string basename(std::string const& path);
std::string extension(std::string const& path);
std::string parent_path(std::string const& path);
bool check_existence(std::string const& path);
void remove_if_exists(std::string const& path);

//todo move to cpp and reduce code duplication!!!
/**
 * Checks if file exists.
 * Analogs: http://www.techbytes.ca/techbyte103.html , http://www.gamedev.net/topic/211918-determining-if-a-file-exists-c/
 */
inline bool FileExists(std::string filename) {
    struct stat st_buf;
    return stat(filename.c_str(), &st_buf) == 0 && S_ISREG(st_buf.st_mode);
}

/**
 * Exit(1) if file doesn't exists, writes FATAL log message.
 */
inline void CheckFileExistenceFATAL(std::string filename) {
  VERIFY_MSG(FileExists(filename), "File " << filename << " doesn't exist or can't be read!\n");
}


// doesn't support symlinks
std::string resolve(std::string const& path);
std::string make_relative_path(std::string p, std::string base = current_dir());
}
