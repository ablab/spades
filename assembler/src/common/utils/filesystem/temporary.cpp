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
#include <fcntl.h>

namespace fs {
namespace impl {
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
        : parent_(parent), fd_(-1), released_(false) {
    std::string tprefix = prefix + "_XXXXXX";
    if (parent)
        tprefix = fs::append_path(parent->dir(), tprefix);
    char *tempprefix = strdup(tprefix.c_str());
    VERIFY_MSG(-1 != (fd_ = ::mkstemp(tempprefix)), "Cannot create temporary file");
    file_ = tempprefix;
    free(tempprefix);
    TRACE("Creating " << file_);
}

TmpFileImpl::TmpFileImpl(nullptr_t, const std::string &file, TmpDir parent)
        : file_(file), parent_(parent), fd_(-1), released_(false) {
    fd_ = ::open(file_.c_str(), O_CREAT | O_RDWR, 0600);
    VERIFY_MSG(-1 != fd_, "Cannot open file");
    TRACE("Acquiring " << file_);
}

TmpFileImpl::~TmpFileImpl() {
    close();
    if (!released_) {
        TRACE("Removing " << file_);
        ::unlink(file_.c_str());
    }
}

void TmpFileImpl::close() {
    ::close(fd_);
}

const std::string &TmpFileImpl::release() {
    VERIFY_MSG(!released_.exchange(true), "Temp file is already released");
    return file_;
}

DependentTmpFile TmpFileImpl::CreateDep(const std::string &suffix) {
    return new DependentTmpFileImpl(suffix, this);
}

DependentTmpFileImpl::DependentTmpFileImpl(const std::string &suffix, TmpFile parent)
        : parent_(parent), file_(parent_->file() + "." + suffix), released_(false) {
    TRACE("Dependent: " << file_);
}

DependentTmpFileImpl::~DependentTmpFileImpl() {
    if (!released_) {
        TRACE("Removing " << file_);
        ::unlink(file_.c_str());
    }
}

const std::string &DependentTmpFileImpl::release() {
    VERIFY_MSG(!released_.exchange(true), "Temp file is already released");
    return file_;
}

}  // namespace impl
}  // namespace fs
