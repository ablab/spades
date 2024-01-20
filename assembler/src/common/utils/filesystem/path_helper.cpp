//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "path_helper.hpp"

#include "utils/logger/logger.hpp"
#include "utils/verify.hpp"

#include <sys/stat.h>
#include <sys/types.h>

#include <dirent.h>
#include <unistd.h>

#include <cstring>

namespace fs {

std::filesystem::path make_temp_dir(std::filesystem::path const& prefix,
                          std::string const& suffix) {
    std::filesystem::path name = prefix / (suffix + "_XXXXXX");
    char* actual;
    if ((actual = ::mkdtemp(strcpy(new char[name.native().length() + 1], name.c_str())))
            == NULL)
        throw std::runtime_error("Cannot create temporary dir " + name.native());

    std::filesystem::path result(actual);
    if (result == name)
        throw std::runtime_error("Cannot create temporary dir " + name.native());

    delete[] actual;

    return result;
}

std::filesystem::path make_full_path(std::filesystem::path const& path) {
    return std::filesystem::current_path() / path;
}

//TODO do we need to screen anything but whitespaces?
std::filesystem::path screen_whitespaces(std::filesystem::path const &path) {
    std::string res = "";
    for (size_t i = 0; i < path.native().size(); i++) {
        if ((i == 0) || (path.native()[i] != ' ') || (path.native()[i - 1] == '\\')) {
            res += path.native()[i];
        } else {
            res +='\\';
            res +=' ';
        }
    }
    return res;
}


std::filesystem::path resolve(const std::filesystem::path &path) {
    std::filesystem::path absPath = absolute(path);
    std::filesystem::path::iterator it = absPath.begin();
    std::filesystem::path result = *it++;

    for (; exists(result / *it) && it != absPath.end(); ++it) {
        result /= *it;
    }
    result = canonical(result);

    for (; it != absPath.end(); ++it) {
        if (*it == "..") {
            result = result.parent_path();
        }
        else if (*it != ".") {
            result /= *it;
        }
    }
    return result.make_preferred();
}
}
