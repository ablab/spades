//***************************************************************************
//* Copyright (c) 2017 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "memory_limit.hpp"

#include "utils/logger/logger.hpp"
#include "utils/verify.hpp"

#if __DARWIN || __DARWIN_UNIX03
# include <mach/task.h>
# include <mach/mach.h>
#else
# include <sys/resource.h>
#endif

#include <sys/time.h>
#include <sys/resource.h>

#include "config.hpp"

#ifdef SPADES_USE_JEMALLOC
# include <jemalloc/jemalloc.h>
#endif

#include <common/utils/logger/logger.hpp>

namespace utils {

void limit_memory(size_t limit) {
    rlimit rl;
    if (sizeof(rlim_t) < 8) {
        FATAL_ERROR("Can't limit virtual memory because of 32-bit system");
        return;
    }

    int res = getrlimit(RLIMIT_AS, &rl);
    if (res != 0)
        FATAL_ERROR("getrlimit(2) call failed, errno = " << errno);

    // We cannot go beyond hard limit and we might not have enough privileges to
    // increase the hard limit
    rl.rlim_cur = std::min<size_t>(limit, rl.rlim_max);
    res = setrlimit(RLIMIT_AS, &rl);
    if (res != 0)
        FATAL_ERROR("setrlimit(2) call failed, errno = " << errno);
    INFO("Memory limit set to " << (1.0 * (double) rl.rlim_cur / 1024 / 1024 / 1024) << " Gb");
}

size_t get_memory_limit() {
    rlimit rl;
    int res = getrlimit(RLIMIT_AS, &rl);
    if (res != 0)
        FATAL_ERROR("getrlimit(2) call failed, errno = " << errno);

    return rl.rlim_cur;
}

#if __DARWIN || __DARWIN_UNIX03
size_t get_max_rss() {
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

  if (KERN_SUCCESS !=
      task_info(mach_task_self(),
                TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count))
    return -1U;

  return t_info.resident_size / 1024;
}
#else

size_t get_max_rss() {
    rusage ru;
    getrusage(RUSAGE_SELF, &ru);

    return ru.ru_maxrss;
}

#endif

size_t get_used_memory() {
#ifdef SPADES_USE_JEMALLOC
    const size_t *cmem = 0;
    size_t clen = sizeof(cmem);

    je_mallctl("stats.cactive", &cmem, &clen, NULL, 0);
    return *cmem;
#else
    return get_max_rss();
#endif
}

size_t get_free_memory() {
    return get_memory_limit() - get_used_memory();
}

}
