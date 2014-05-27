//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
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
#include "logger/logger.hpp"
#include "verify.hpp"

namespace path {
//todo review and make names consistent!

typedef std::vector<std::string> files_t;

bool make_dir(std::string const& folder);
std::string make_temp_dir(std::string const& prefix, std::string const& suffix);
void remove_dir(std::string const& folder);
bool is_regular_file(std::string const& path);
std::string append_path(std::string const& prefix, std::string const& suffix);
std::string current_dir();

//todo why non-cons argument?!
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
  if(!FileExists(filename))
      FATAL_ERROR("File " << filename << " doesn't exist or can't be read!");
}

inline void make_dirs(const std::string& path) {
  VERIFY(!path.empty());

  size_t slash_pos = 0;
  while ((slash_pos = path.find_first_of('/', slash_pos + 1)) != std::string::npos) {
    make_dir(path.substr(0, slash_pos));
  }
  if (path[path.size() - 1] != '/') {
    make_dir(path);
  }
}

// doesn't support symlinks
std::string resolve(std::string const& path);
std::string make_relative_path(std::string p, std::string base = current_dir());
}
