//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "utils/verify.hpp"

namespace utils {

inline rlim_t limit_file(size_t limit) {
  struct rlimit rl;

  int res = getrlimit(RLIMIT_NOFILE, &rl);
  VERIFY_MSG(res == 0,
             "getrlimit(2) call failed, errno = " << errno);

  // We cannot go beyond hard limit and we might not have enough privileges to
  // increase the hard limit
  limit = std::max<size_t>(limit, rl.rlim_cur);
  rl.rlim_cur = std::min<size_t>(limit, rl.rlim_max);
  res = setrlimit(RLIMIT_NOFILE, &rl);
  VERIFY_MSG(res == 0,
             "setrlimit(2) call failed, errno = " << errno);
  INFO("Open file limit set to " << rl.rlim_cur);

  return rl.rlim_cur;
}

}
