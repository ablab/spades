//===-- Debug.cpp - An easy way to add debug output to your code ----------===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// This file implements a handy way of adding debugging information to your
// code, without it being enabled all of the time, and without having to add
// command line options to enable it.
//
// In particular, just wrap your code with the LLVM_DEBUG() macro, and it will
// be enabled automatically if you specify '-debug' on the command-line.
// Alternatively, you can also use the SET_DEBUG_TYPE("foo") macro to specify
// that your debug code belongs to class "foo".  Then, on the command line, you
// can specify '-debug-only=foo' to enable JUST the debug information for the
// foo class.
//
// When compiling without assertions, the -debug-* options and all code in
// LLVM_DEBUG() statements disappears, so it does not affect the runtime of the
// code.
//
//===----------------------------------------------------------------------===//

#include "llvm/Support/Debug.h"
#include "llvm/Support/ManagedStatic.h"
#include "llvm/Support/Signals.h"
#include "llvm/Support/raw_ostream.h"

#include <vector>
#include <string>

#undef isCurrentDebugType
#undef setCurrentDebugType
#undef setCurrentDebugTypes

using namespace llvm;

// Even though LLVM might be built with NDEBUG, define symbols that the code
// built without NDEBUG can depend on via the llvm/Support/Debug.h header.
namespace llvm {
/// Exported boolean set by the -debug option.
bool DebugFlag = false;

static ManagedStatic<std::vector<std::string>> CurrentDebugType;

/// Return true if the specified string is the debug type
/// specified on the command line, or if none was specified on the command line
/// with the -debug-only=X option.
bool isCurrentDebugType(const char *DebugType) {
  if (CurrentDebugType->empty())
    return true;
  // See if DebugType is in list. Note: do not use find() as that forces us to
  // unnecessarily create an std::string instance.
  for (auto &d : *CurrentDebugType) {
    if (d == DebugType)
      return true;
  }
  return false;
}

/// Set the current debug type, as if the -debug-only=X
/// option were specified.  Note that DebugFlag also needs to be set to true for
/// debug output to be produced.
///
void setCurrentDebugTypes(const char **Types, unsigned Count);

void setCurrentDebugType(const char *Type) {
  setCurrentDebugTypes(&Type, 1);
}

void setCurrentDebugTypes(const char **Types, unsigned Count) {
  CurrentDebugType->clear();
  for (size_t T = 0; T < Count; ++T)
    CurrentDebugType->push_back(Types[T]);
}
} // namespace llvm

// Avoid "has no symbols" warning.
namespace llvm {
  /// dbgs - Return errs().
  raw_ostream &dbgs() {
    return errs();
  }
}

/// EnableDebugBuffering - Turn on signal handler installation.
///
bool llvm::EnableDebugBuffering = false;
