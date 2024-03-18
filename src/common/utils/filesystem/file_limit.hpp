//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/verify.hpp"
#include "utils/logger/logger.hpp"

#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <limits.h>

#include <cstring>

namespace utils {

inline rlim_t limit_file(size_t limit) {
  struct rlimit rl;

  int res = getrlimit(RLIMIT_NOFILE, &rl);
  CHECK_FATAL_ERROR(res == 0,
             "getrlimit(2) call failed, errno = " << errno);

  // We cannot go beyond hard limit and we might not have enough privileges to
  // increase the hard limit
  limit = std::max<size_t>(limit, rl.rlim_cur);
  rl.rlim_cur = std::min<size_t>(limit, rl.rlim_max);
  // If OPEN_MAX is defined, then limit by it as well
#ifdef OPEN_MAX
  rl.rlim_cur = std::min<size_t>(limit, OPEN_MAX);
#endif
  res = setrlimit(RLIMIT_NOFILE, &rl);
  if (res != 0) {
      WARN("Failed to set file limit to " << rl.rlim_cur << ", setrlimit(2) call failed, errno = "
           << errno << " (" << strerror(errno) << ")");
  } else {
      INFO("Open file limit set to " << rl.rlim_cur);
  }
  
  return rl.rlim_cur;
}

}
