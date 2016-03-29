//===- llvm/Support/Unix/Unix.h - Common Unix Include File -------*- C++ -*-===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
//
// This file defines things specific to Unix implementations.
//
//===----------------------------------------------------------------------===//

#ifndef LLVM_LIB_SUPPORT_UNIX_UNIX_H
#define LLVM_LIB_SUPPORT_UNIX_UNIX_H

//===----------------------------------------------------------------------===//
//=== WARNING: Implementation here must contain only generic UNIX code that
//===          is guaranteed to work on all UNIX variants.
//===----------------------------------------------------------------------===//

#include "llvm/Support/Errno.h"
#include <algorithm>
#include <assert.h>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>

#include <unistd.h>
#include <sys/param.h>
#include <sys/time.h>
#include <time.h>

#include <dlfcn.h>

/// This function builds an error message into \p ErrMsg using the \p prefix
/// string and the Unix error number given by \p errnum. If errnum is -1, the
/// default then the value of errno is used.
/// @brief Make an error message
///
/// If the error number can be converted to a string, it will be
/// separated from prefix by ": ".
static inline bool MakeErrMsg(
  std::string* ErrMsg, const std::string& prefix, int errnum = -1) {
  if (!ErrMsg)
    return true;
  if (errnum == -1)
    errnum = errno;
  *ErrMsg = prefix + ": " + llvm::sys::StrError(errnum);
  return true;
}

#endif
