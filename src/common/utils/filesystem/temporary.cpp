//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2017-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "temporary.hpp"

#include "io/binary/binary.hpp"
#include "io/binary/types/fs.hpp"
#include "path_helper.hpp"
#include "utils/verify.hpp"
#include "utils/logger/logger.hpp"

#include <cerrno>
#include <cstring>
#include <unistd.h>
#include <fcntl.h>

namespace fs {
namespace impl {
TmpDirImpl::TmpDirImpl(const std::filesystem::path &prefix, const std::string &suffix)
        : dir_(fs::make_temp_dir(prefix, suffix)), released_(false) {
    TRACE("Creating " << dir_);
}

TmpDirImpl::TmpDirImpl(std::nullptr_t, const std::filesystem::path &dir)
        : dir_(dir), released_(false) {
    TRACE("Acquiring " << dir_);
}

TmpDirImpl::~TmpDirImpl() {
    if (released_)
        return;

    TRACE("Removing " << dir_);
    remove_all(dir_);
}

const std::filesystem::path &TmpDirImpl::release() {
    bool already_released = released_.exchange(true);
    CHECK_FATAL_ERROR(!already_released, "Temp dir is already released");
    return dir_;
}

TmpFile TmpDirImpl::tmp_file(const std::string &prefix) {
    return new TmpFileImpl(prefix, this);
}

TmpFileImpl::TmpFileImpl(const std::string &prefix, TmpDir parent)
        : parent_(parent), fd_(-1), released_(false) {
    std::filesystem::path tprefix = prefix + "_XXXXXX";
    if (parent)
        tprefix = parent->dir() / tprefix;
    char *tempprefix = strdup(tprefix.c_str());
    CHECK_FATAL_ERROR(-1 != (fd_ = ::mkstemp(tempprefix)), "Cannot create temporary file");
    file_ = tempprefix;
    free(tempprefix);
    TRACE("Creating " << file_);
}

TmpFileImpl::TmpFileImpl(std::nullptr_t, const std::filesystem::path &file, TmpDir parent)
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

const std::filesystem::path &TmpFileImpl::release() {
    CHECK_FATAL_ERROR(!released_.exchange(true), "Temp file is already released");
    return file_;
}

DependentTmpFile TmpFileImpl::CreateDep(const std::string &suffix) {
    return new DependentTmpFileImpl(suffix, this);
}

DependentTmpFileImpl::DependentTmpFileImpl(const std::string &suffix, TmpFile parent)
        : parent_(parent), file_(parent_->file().native() + "." + suffix), released_(false) {
    TRACE("Dependent: " << file_);
}

DependentTmpFileImpl::DependentTmpFileImpl(std::nullptr_t, const std::string &file, TmpFile parent)
        : parent_(parent), file_(file), released_(false) {
    TRACE("Dependent: " << file_);
}

DependentTmpFileImpl::~DependentTmpFileImpl() {
    if (released_) {
        return;
    }
    TRACE("Removing " << file_);
    ::unlink(file_.c_str());
}

const std::filesystem::path &DependentTmpFileImpl::release() {
    bool already_released = released_.exchange(true);
    CHECK_FATAL_ERROR(!already_released, "Temp file is already released");
    return file_;
}

}  // namespace impl
}  // namespace fs

namespace io {
namespace binary {
namespace impl {
void Serializer<fs::TmpDir>::Read(std::istream &is, fs::TmpDir &dir) {
    bool initialized = false;
    io::binary::BinRead(is, initialized);
    if (initialized) {
        std::string str;
        io::binary::BinRead(is, str);
        dir = fs::tmp::acquire_temp_dir(str);
        dir->release();
    } else
        dir.reset();
}

void Serializer<fs::TmpDir>::Write(std::ostream &os, const fs::TmpDir &dir) {
    io::binary::BinWrite(os, (bool)dir);
    if (dir)
        io::binary::BinWrite(os, dir->dir());
}

void Serializer<fs::TmpFile>::Read(std::istream &is, fs::TmpFile &f) {
    bool initialized = false;
    io::binary::BinRead(is, initialized);
    if (initialized) {
        std::string str;
        io::binary::BinRead(is, str);
        f = fs::tmp::acquire_temp_file(str);
        f->release();
    } else
        f.reset();
}

void Serializer<fs::TmpFile>::Write(std::ostream &os, const fs::TmpFile &f) {
    io::binary::BinWrite(os, (bool)f);
    if (f)
        io::binary::BinWrite(os, f->file());
}

void Serializer<fs::DependentTmpFile>::Read(std::istream &is, fs::DependentTmpFile &f) {
    bool initialized = false;
    io::binary::BinRead(is, initialized);
    if (initialized) {
        std::string str;
        io::binary::BinRead(is, str);
        f = fs::tmp::acquire_dependent_temp_file(str);
        f->release();
    } else
        f.reset();
}

void Serializer<fs::DependentTmpFile>::Write(std::ostream &os, const fs::DependentTmpFile &f) {
    io::binary::BinWrite(os, (bool)f);
    if (f)
        io::binary::BinWrite(os, f->file());
}
}
}
}
