/**
 * @file timer.h
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * Use CPU instructions and inline assembly to query tics and/or wall clock.
 */
#ifndef __TIMER_H__
#define __TIMER_H__

#ifdef __cplusplus
extern "C" {
#endif

/*#define FORCE_SYS_TIME 1*/

#if (defined(__i386__) || defined(__x86_64__) || defined(__powerpc__)) && !defined(_CRAYC) && !defined(FORCE_SYS_TIME)
#   define HAVE_RDTSC 1
#   if defined(__i386__)
static __inline__ unsigned long long rdtsc(void)
{
    unsigned long long int x;
    __asm__ volatile(".byte 0x0f, 0x31" : "=A"(x));
    return x;
}
#   elif defined(__x86_64__)
static __inline__ unsigned long long rdtsc(void)
{
    unsigned hi, lo;
    __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
    return ((unsigned long long)lo) | (((unsigned long long)hi) << 32);
}
#   elif defined(__powerpc__)
static __inline__ unsigned long long rdtsc(void)
{
    unsigned long long int result = 0;
    unsigned long int upper, lower, tmp;
    __asm__ volatile(
        "0:                  \n"
        "\tmftbu   %0        \n"
        "\tmftb    %1        \n"
        "\tmftbu   %2        \n"
        "\tcmpw    %2,%0     \n"
        "\tbne     0b        \n"
        : "=r"(upper), "=r"(lower), "=r"(tmp)
    );
    result = upper;
    result = result << 32;
    result = result | lower;

    return (result);
}
#   endif
#elif defined(HAVE_SYS_TIME_H) || defined(FORCE_SYS_TIME)
#   undef HAVE_SYS_TIME_H
#   define HAVE_SYS_TIME_H 1
#   include <sys/time.h>
#elif defined(HAVE_WINDOWS_H)
#   undef HAVE_WINDOWS_H
#   define HAVE_WINDOWS_H 1
#   include <windows.h>
static LARGE_INTEGER frequency;
#endif

static unsigned long long timer_start()
{
#if defined(HAVE_RDTSC)
    return rdtsc();
#elif defined(HAVE_SYS_TIME_H)
    struct timeval timer;
    (void)gettimeofday(&timer, NULL);
    return timer.tv_sec * 1000000 + timer.tv_usec;
#elif defined(HAVE_WINDOWS_H)
    LARGE_INTEGER timer;
    QueryPerformanceCounter(&timer);
    return timer.QuadPart * 1000 / frequency.QuadPart;
#else
    return 0;
#endif
}

/*
static unsigned long long timer_end(unsigned long long begin)
{
    return timer_start() - begin;
}
*/

/*
static void timer_init()
{
#if defined(HAVE_RDTSC)
#elif defined(HAVE_SYS_TIME_H)
#elif defined(HAVE_WINDOWS_H)
    QueryPerformanceFrequency(&frequency);
#else
#endif
}
*/

/*
static const char *timer_name()
{
#if defined(HAVE_RDTSC)
    return "rdtsc";
#elif defined(HAVE_SYS_TIME_H)
    return "gettimeofday";
#elif defined(HAVE_WINDOWS_H)
    return "windows QueryPerformanceCounter";
#else
    return "no timers";
#endif
}
*/

#ifdef __cplusplus
}
#endif

#endif /* __TIMER_H__ */
