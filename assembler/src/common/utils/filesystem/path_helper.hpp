//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <string>
#include <vector>

namespace fs {
//todo review and make names consistent!

typedef std::vector<std::string> files_t;

bool make_dir(std::string const &folder);

std::string make_temp_dir(std::string const &prefix, std::string const &suffix);

void remove_dir(std::string const &folder);

bool is_regular_file(std::string const &path);

std::string append_path(std::string const &prefix, std::string const &suffix);

std::string current_dir();

std::string make_full_path(std::string const &path);

std::string filename(std::string const &path);

std::string basename(std::string const &path);

std::string extension(std::string const &path);

std::string parent_path(std::string const &path);

bool check_existence(std::string const &path);

void remove_if_exists(std::string const &path);

std::string screen_whitespaces(std::string const &path);

/**
* Checks if file exists.
* Analogs: http://www.techbytes.ca/techbyte103.html , http://www.gamedev.net/topic/211918-determining-if-a-file-exists-c/
*/
bool FileExists(std::string const &filename);

/**
* Exit(1) if file doesn't exists, writes FATAL log message.
*/
void CheckFileExistenceFATAL(std::string const &filename);

void make_dirs(std::string const &path);

// doesn't support symlinks
std::string resolve(std::string const &path);

std::string make_relative_path(std::string p, std::string base = current_dir());

std::string MakeLaunchTimeDirName();

class TmpFolderFixture
{
    std::string tmp_folder_;

public:
    TmpFolderFixture(std::string tmp_folder = "tmp") :
            tmp_folder_(tmp_folder)
    {
        fs::make_dirs(tmp_folder_);
    }

    ~TmpFolderFixture()
    {
        fs::remove_dir(tmp_folder_);
    }
};

}
