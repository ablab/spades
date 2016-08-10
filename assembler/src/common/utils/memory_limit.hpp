//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#if __DARWIN || __DARWIN_UNIX03
#include <mach/task.h>
#include <mach/mach.h>
#else

#include <sys/resource.h>

#endif

#include <sys/time.h>
#include <sys/resource.h>

#include "config.hpp"

#ifdef SPADES_USE_JEMALLOC

# include <jemalloc/jemalloc.h>

#endif

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
    INFO("Memory limit set to " << (1.0 * (double) rl.rlim_cur / 1024 / 1024 / 1024) << " Gb");
}

inline size_t get_memory_limit() {
    rlimit rl;
    int res = getrlimit(RLIMIT_AS, &rl);
    VERIFY_MSG(res == 0,
               "getrlimit(2) call failed, errno = " << errno);

    return rl.rlim_cur;
}

#if __DARWIN || __DARWIN_UNIX03
inline size_t get_max_rss() {
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

  if (KERN_SUCCESS !=
      task_info(mach_task_self(),
                TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count))
    return -1U;

  return t_info.resident_size / 1024;
}
#else

inline size_t get_max_rss() {
    rusage ru;
    getrusage(RUSAGE_SELF, &ru);

    return ru.ru_maxrss;
}

#endif

inline size_t get_used_memory() {
#ifdef SPADES_USE_JEMALLOC
    const size_t *cmem = 0;
    size_t clen = sizeof(cmem);

    je_mallctl("stats.cactive", &cmem, &clen, NULL, 0);
    return *cmem;
#else
    get_max_rss();
#endif
}


inline size_t get_free_memory() {
    return get_memory_limit() - get_used_memory();
}
