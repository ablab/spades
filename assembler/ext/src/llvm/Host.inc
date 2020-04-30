//===- llvm/Support/Unix/Host.inc -------------------------------*- C++ -*-===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// This file implements the UNIX Host support.
//
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
//=== WARNING: Implementation here must contain only generic UNIX code that
//===          is guaranteed to work on *all* UNIX variants.
//===----------------------------------------------------------------------===//

#include "Unix.h"
#include "llvm/ADT/StringRef.h"
#include "llvm/Config.h"
#include <cctype>
#include <string>
#include <sys/utsname.h>

using namespace llvm;

static std::string getOSVersion() {
  struct utsname info;

  if (uname(&info))
    return "";

  return info.release;
}