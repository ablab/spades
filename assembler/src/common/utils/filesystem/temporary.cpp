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
        : dir_(fs::make_temp_dir(prefix, suffix)), released_(false) {
    TRACE("Creating " << dir_);
}

TmpDirImpl::TmpDirImpl(std::nullptr_t, const std::string &dir)
        : dir_(dir), released_(false) {
    TRACE("Acquiring " << dir_);
}

TmpDirImpl::~TmpDirImpl() {
    if (released_) {
        return;
    }
    TRACE("Removing " << dir_);
    fs::remove_dir(dir_);
}

const std::string &TmpDirImpl::release() {
    bool already_released = released_.exchange(true);
    CHECK_FATAL_ERROR(!already_released, "Temp dir is already released");
    return dir_;
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
    CHECK_FATAL_ERROR(-1 != (fd_ = ::mkstemp(tempprefix)), "Cannot create temporary file");
    file_ = tempprefix;
    free(tempprefix);
    TRACE("Creating " << file_);
}

TmpFileImpl::TmpFileImpl(std::nullptr_t, const std::string &file, TmpDir parent)
        : file_(file), parent_(parent), fd_(-1), released_(false) {
    fd_ = ::open(file_.c_str(), O_CREAT | O_RDWR, 0600);
    CHECK_FATAL_ERROR(-1 != fd_, "Cannot open file");
    TRACE("Acquiring " << file_);
}

TmpFileImpl::~TmpFileImpl() {
    close();
    if (released_) {
        return;
    }
    TRACE("Removing " << file_);
    ::unlink(file_.c_str());
}

void TmpFileImpl::close() {
    ::close(fd_);
}

const std::string &TmpFileImpl::release() {
    CHECK_FATAL_ERROR(!released_.exchange(true), "Temp file is already released");
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
    if (released_) {
        return;
    }
    TRACE("Removing " << file_);
    ::unlink(file_.c_str());
}

const std::string &DependentTmpFileImpl::release() {
    bool already_released = released_.exchange(true);
    CHECK_FATAL_ERROR(!already_released, "Temp file is already released");
    return file_;
}

}  // namespace impl
}  // namespace fs
