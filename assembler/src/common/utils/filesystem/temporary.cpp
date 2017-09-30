//***************************************************************************
//* Copyright (c) 2017 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "temporary.hpp"

#include "path_helper.hpp"
#include "utils/verify.hpp"
#include "utils/logger/logger.hpp"

#include <cerrno>
#include <cstring>
#include <unistd.h>

namespace fs {
TmpDirImpl::TmpDirImpl(const std::string &prefix, const std::string &suffix)
        : dir_(fs::make_temp_dir(prefix, suffix)) {
    TRACE("Creating " << dir_);
}

TmpDirImpl::~TmpDirImpl() {
    TRACE("Removing " << dir_);
    fs::remove_dir(dir_);
}

TmpFile TmpDirImpl::tmp_file(const std::string &prefix) {
    return new TmpFileImpl(prefix, this);
}

TmpFileImpl::TmpFileImpl(const std::string &prefix, TmpDir parent)
        : parent_(parent), fd_(-1) {
    std::string tprefix = prefix + "_XXXXXX";
    if (parent)
        tprefix = fs::append_path(parent->dir(), tprefix);
    char *tempprefix = strdup(tprefix.c_str());
    VERIFY_MSG(-1 != (fd_ = ::mkstemp(tempprefix)), "Cannot create temporary file");
    file_ = tempprefix;
    free(tempprefix);
    TRACE("Creating " << file_);
}

TmpFileImpl::~TmpFileImpl() {
    TRACE("Removing " << file_);
    close();
    ::unlink(file_.c_str());
}

void TmpFileImpl::close() {
    ::close(fd_);
}

DependentTmpFile TmpFileImpl::CreateDep(const std::string &suffix) {
    return new DependentTmpFileImpl(suffix, this);
}

DependentTmpFileImpl::DependentTmpFileImpl(const std::string &suffix, TmpFile parent)
        : parent_(parent), file_(parent_->file() + "." + suffix) {
    TRACE("Dependent: " << file_);
}

DependentTmpFileImpl::~DependentTmpFileImpl() {
    TRACE("Removing " << file_);
    ::unlink(file_.c_str());
}

}
