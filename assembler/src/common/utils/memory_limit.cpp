//***************************************************************************
//* Copyright (c) 2017 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "memory_limit.hpp"

#include "utils/parallel/openmp_wrapper.h"
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
    return -1ULL;

  return t_info.resident_size / 1024;
}
#else

size_t get_max_rss() {
    rusage ru;
    getrusage(RUSAGE_SELF, &ru);

    return ru.ru_maxrss;
}

#endif

#if defined(SPADES_USE_MIMALLOC)
extern "C" {
    void mi_stats_merge(void);
    void mi_collect(bool);
    size_t mi_stats_total_mem();
};
#endif

size_t get_used_memory() {
#if defined(SPADES_USE_JEMALLOC)
    // Update statistics cached by mallctl
    {
        uint64_t epoch = 1;
        size_t sz = sizeof(epoch);
        if (je_mallctl("epoch", &epoch, &sz, &epoch, sz) != 0)
            FATAL_ERROR("mallctl() call failed, errno = " << errno);
    }

    {
        size_t cmem = 0;
        size_t clen = sizeof(cmem);

        int res = je_mallctl("stats.active", &cmem, &clen, NULL, 0);
        if (res != 0)
            FATAL_ERROR("mallctl() call failed, errno = " << errno);

        return cmem;
    }
#elif defined(SPADES_USE_MIMALLOC)
    // mimalloc implements separate and independent memory pulls for each thread
    // The statistics is also collected per pool. So we essentially need to propagate
    // the stats from per-thread pool into main one
    if (omp_get_thread_num() > 0) {
        mi_stats_merge();
    } else {
        unsigned nthreads = omp_get_max_threads();
#       pragma omp parallel for
        for (unsigned i = 0; i < 2*nthreads; ++i) {
            mi_collect(true); // FIXME: hack-hack-hack
            mi_stats_merge();
        }
    }
    return mi_stats_total_mem();
#else
    return get_max_rss();
#endif
}

size_t get_free_memory() {
    return get_memory_limit() - get_used_memory();
}

}
