//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include <sys/time.h>
#include <sys/resource.h>

inline void limit_memory(size_t limit) {
  rlimit rl;
  if (sizeof(rlim_t) < 8) {
    INFO("Can't limit virtual memory because of 32-bit system");
    return;
  }

  int res = getrlimit(RLIMIT_AS, &rl);
  VERIFY_MSG(res == 0,
             "getrlimit(2) call failed, errno = " << errno);

  // We cannot go beyond hard limit and we might not have enough privileges to
  // increase the hard limit
  rl.rlim_cur = std::min<size_t>(limit, rl.rlim_max);
  res = setrlimit(RLIMIT_AS, &rl);
  VERIFY_MSG(res == 0,
             "setrlimit(2) call failed, errno = " << errno);
  INFO("Memory limit set to " << (1.0 * (double)rl.rlim_cur / 1024 / 1024 / 1024) << " Gb");
}

inline size_t get_memory_limit() {
  rlimit rl;
  int res = getrlimit(RLIMIT_AS, &rl);
  VERIFY_MSG(res == 0,
             "getrlimit(2) call failed, errno = " << errno);

  return rl.rlim_cur;
}
