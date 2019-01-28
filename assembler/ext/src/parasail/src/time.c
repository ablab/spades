/**
 * @file timer_real.h
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * Portable use of nanosecond precision realtime clock.
 */
#include "config.h"

#include "parasail/parasail.h"

#if HAVE_GETSYSTEMTIMEASFILETIME && !defined(INT64_LITERAL_SUFFIX_UNKNOWN)
#include <windows.h>
#elif HAVE_CLOCK_GET_TIME
#include <mach/clock.h>
#include <mach/mach.h>
#elif HAVE_CLOCK_GETTIME_MONOTONIC || HAVE_CLOCK_GETTIME_REALTIME
#define _POSIX_C_SOURCE 199309L
#include <time.h>
#endif

#include <assert.h>

double parasail_time(void)
{
#if HAVE_GETSYSTEMTIMEASFILETIME && !defined(INT64_LITERAL_SUFFIX_UNKNOWN)
    __int64 wintime;
    double sec;
    double nsec;
    GetSystemTimeAsFileTime((FILETIME*)&wintime);
#if defined(INT64_LITERAL_SUFFIX_I64)
    wintime -=116444736000000000i64;   /*1jan1601 to 1jan1970*/
    sec  = wintime / 10000000i64;      /*seconds*/
    nsec = wintime % 10000000i64 *100; /*nano-seconds*/
#elif defined(INT64_LITERAL_SUFFIX_LL)
    wintime -=116444736000000000LL;   /*1jan1601 to 1jan1970*/
    sec  = wintime / 10000000LL;      /*seconds*/
    nsec = wintime % 10000000LL *100; /*nano-seconds*/
#endif
    return sec + nsec/1000000000.0;
#elif HAVE_CLOCK_GET_TIME
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    return (double)(mts.tv_sec) + (double)(mts.tv_nsec)/1000000000.0;
#elif HAVE_CLOCK_GETTIME_MONOTONIC
    struct timespec ts;
    /* Works on FreeBSD */
    long retval = clock_gettime(CLOCK_MONOTONIC, &ts);
    assert(0 == retval);
    return (double)(ts.tv_sec) + (double)(ts.tv_nsec)/1000000000.0;
#elif HAVE_CLOCK_GETTIME_REALTIME
    struct timespec ts;
    /* Works on Linux */
    long retval = clock_gettime(CLOCK_REALTIME, &ts);
    assert(0 == retval);
    return (double)(ts.tv_sec) + (double)(ts.tv_nsec)/1000000000.0;
#endif
    return 0.0; /* avoid warnings */
}

