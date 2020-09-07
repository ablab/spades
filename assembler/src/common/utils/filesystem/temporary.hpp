//***************************************************************************
//* Copyright (c) 2017 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/logger/decl_logger.hpp"

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

class non_copy_move_assign_able {
  protected:
    non_copy_move_assign_able() = default;
    ~non_copy_move_assign_able() = default;

  private:
    non_copy_move_assign_able(const non_copy_move_assign_able&) = delete;
    non_copy_move_assign_able(non_copy_move_assign_able&&) = delete;
    non_copy_move_assign_able &operator=(const non_copy_move_assign_able&) = delete;
    non_copy_move_assign_able &operator=(non_copy_move_assign_able&&) = delete;
};

class TmpDirImpl : public llvm::ThreadSafeRefCountedBase<TmpDirImpl>, non_copy_move_assign_able {
  public:
    // Create new tmp dir
    TmpDirImpl(const std::string &prefix, const std::string &suffix);
    // Acquire existing tmp dir
    TmpDirImpl(std::nullptr_t, const std::string &dir);
    ~TmpDirImpl();

    const std::string &dir() const { return dir_; }
    operator std::string() const { return dir_; }
    TmpFile tmp_file(const std::string &prefix = "tmp");
    const std::string &release();

  private:
    std::string dir_;
    std::atomic<bool> released_;

    DECL_LOGGER("Temporary");
};

class TmpFileImpl : public llvm::ThreadSafeRefCountedBase<TmpFileImpl>, non_copy_move_assign_able {
  public:
    // Create new tmp file
    TmpFileImpl(const std::string &prefix = "tmp", TmpDir parent = nullptr);
    // Acquire existing file
    TmpFileImpl(std::nullptr_t, const std::string &file, TmpDir parent = nullptr);
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

    DECL_LOGGER("Temporary");
};

class DependentTmpFileImpl : public llvm::ThreadSafeRefCountedBase<DependentTmpFileImpl>, non_copy_move_assign_able {
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

    DECL_LOGGER("Temporary");
};

inline TmpDir make_temp_dir(const std::string &prefix, const std::string &suffix) {
    return new TmpDirImpl(prefix, suffix);
}

inline TmpFile make_temp_file(const std::string &prefix = "tmp", TmpDir parent = nullptr) {
    return new TmpFileImpl(prefix, parent);
}

inline TmpDir acquire_temp_dir(const std::string &dir) {
    return new TmpDirImpl(nullptr, dir);
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
using impl::acquire_temp_dir;
using impl::make_temp_dir;
using impl::make_temp_file;
}  // namespace tmp
}  // namespace fs
