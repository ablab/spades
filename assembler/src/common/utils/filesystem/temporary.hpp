//***************************************************************************
//* Copyright (c) 2017 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <llvm/ADT/IntrusiveRefCntPtr.h>
#include <string>
#include <atomic>

namespace fs {
namespace impl {
class TmpDirImpl;
class TmpFileImpl;
class DependentTmpFileImpl;

typedef llvm::IntrusiveRefCntPtr<TmpDirImpl> TmpDir;
typedef llvm::IntrusiveRefCntPtr<TmpFileImpl> TmpFile;
typedef llvm::IntrusiveRefCntPtr<DependentTmpFileImpl> DependentTmpFile;

class TmpDirImpl : public llvm::ThreadSafeRefCountedBase<TmpDirImpl> {
  public:
    TmpDirImpl(const std::string &prefix, const std::string &suffix);
    ~TmpDirImpl();

    const std::string &dir() const { return dir_; }
    operator std::string() const { return dir_; }
    TmpFile tmp_file(const std::string &prefix = "tmp");

  private:
    std::string dir_;
};

class TmpFileImpl : public llvm::ThreadSafeRefCountedBase<TmpFileImpl> {
  public:
    // Create new tmp file
    TmpFileImpl(const std::string &prefix = "tmp", TmpDir parent = nullptr);
    // Acquire existing file
    TmpFileImpl(nullptr_t, const std::string &file, TmpDir parent = nullptr);
    ~TmpFileImpl();

    const std::string &file() const { return file_; }
    operator std::string() const { return file_; }
    const std::string &dir() const {
        static std::string noparent("");
        return (parent_ ? parent_->dir() : noparent);
    }
    int fd() const { return fd_; }
    void close();
    DependentTmpFile CreateDep(const std::string &suffix);
    const std::string &release();

  private:
    std::string file_;
    TmpDir parent_;
    int fd_;
    std::atomic<bool> released_;
};

class DependentTmpFileImpl : public llvm::ThreadSafeRefCountedBase<DependentTmpFileImpl> {
  public:
    DependentTmpFileImpl(const std::string &suffix, TmpFile parent);
    ~DependentTmpFileImpl();

    const std::string &file() const { return file_; }
    const std::string &dir() const { return parent_->dir(); }
    operator std::string() const { return file_; }
    const std::string &release();

  private:
    TmpFile parent_;
    std::string file_;
    std::atomic<bool> released_;
};

inline TmpDir make_temp_dir(const std::string &prefix, const std::string &suffix) {
    return new TmpDirImpl(prefix, suffix);
}

inline TmpFile make_temp_file(const std::string &prefix = "tmp", TmpDir parent = nullptr) {
    return new TmpFileImpl(prefix, parent);
}

inline TmpFile acquire_temp_file(const std::string &file, TmpDir parent = nullptr) {
    return new TmpFileImpl(nullptr, file, parent);
}

}  // namespace impl

using impl::DependentTmpFile;
using impl::TmpDir;
using impl::TmpFile;

namespace tmp {
using impl::acquire_temp_file;
using impl::make_temp_dir;
using impl::make_temp_file;
}  // namespace tmp
}  // namespace fs
