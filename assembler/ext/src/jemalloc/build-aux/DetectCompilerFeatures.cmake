# This file configures features based on compiler and library availability
include(CheckLibraryExists)
include(CheckSymbolExists)
include(CheckFunctionExists)

include(CheckCCompilerFlag)
include(CheckCSourceCompiles)
include(CheckCSourceRuns)

include(CheckIncludeFiles)

include(TestBigEndian)

check_include_files(alloca.h JEMALLOC_HAS_ALLOCA_H)

check_library_exists(m log "math.h" has_libm)

check_include_files("pthread.h" has_pthread_h)
if(has_pthread_h)
    set(JEMALLOC_HAVE_PTHREAD 1)
    check_library_exists(pthread pthread_create "pthread.h" has_libpthread)
    if(NOT has_libpthread)
        check_library_exists(c pthread_create "pthread.h" has_libc_pthread)
        if(NOT has_libc_pthread)
            message(WARNING "libpthread is missing")
            unset(JEMALLOC_HAVE_PTHREAD)
        endif()
    endif()

    if(JEMALLOC_HAVE_PTHREAD)
        list(APPEND wrap_syms pthread_create)
    endif()
else()
    message(WARNING "pthread.h is missing")
endif()

# Whether malloc_usable_size definition can use const argument
check_include_files(malloc.h has_malloc_h)
set(JEMALLOC_USABLE_SIZE_CONST " ")
if(has_malloc_h)
    set(CMAKE_REQUIRED_FLAGS "-Werror -Wall -Wextra -std=gnu11")
    check_c_source_compiles("
    #include <malloc.h>
    #include <stddef.h>
    size_t malloc_usable_size(const void *ptr);
    int main(void) {
        return 0;
    }" malloc_usable_size_const)

    if (malloc_usable_size_const)
        set(JEMALLOC_USABLE_SIZE_CONST const)
    endif()

    set(CMAKE_REQUIRED_FLAGS)
endif()

test_big_endian(JEMALLOC_BIG_ENDIAN)

set(CMAKE_REQUIRED_FLAGS "-Werror -Wall -Wextra -std=gnu11")
check_c_source_compiles("
__attribute__((unused)) void foo(void){}
int main() {
    foo();
    return 0;
}" JEMALLOC_HAVE_ATTR)
# configure.ac has a section I don't see the relevance of here:
#   if test "x${GCC}" = "xyes" -a "x${abi}" = "xelf"; then
#    JE_CFLAGS_ADD([-fvisibility=hidden])
#    JE_CXXFLAGS_ADD([-fvisibility=hidden])
#  fi
# If those conditions are still necessary, re-implement here too.
check_c_source_compiles("
int main(void) {
    static __thread int __attribute__((tls_model(\"initial-exec\"), unused)) f;
    f = 0;
    return 0;
}" HAVE_TLS_MODEL)
if(HAVE_TLS_MODEL)
    set(JEMALLOC_TLS_MODEL "__attribute__((tls_model(\"initial-exec\")))")
else()
    set(JEMALLOC_TLS_MODEL " ")
endif()

check_c_source_compiles("
#include <stdlib.h>
void *foo(size_t size) __attribute__((alloc_size(1)));
int main(void) {
    return 0;
}" JEMALLOC_HAVE_ATTR_ALLOC_SIZE)

check_c_source_compiles("
#include <stdlib.h>
void *foo(const char *format, ...) __attribute__((format(gnu_printf, 1, 2)));
int main(void) {
    return 0;
}" JEMALLOC_HAVE_ATTR_FORMAT_GNU_PRINTF)

check_c_source_compiles("
#include <stdlib.h>
void *foo(const char *format, ...) __attribute__((format(printf, 1, 2)));
__attribute__((format(printf, 1, 2))) void *foo(const char *format, ...) {
    return (void *)format;
}
int main() {
    (void)foo(NULL);
    return 0;
}" JEMALLOC_HAVE_ATTR_FORMAT_PRINTF)

check_c_source_compiles("
#include <stddef.h>
int main() {
    void *restrict abc = NULL;
    (void)abc;
    return 0;
}" JEMALLOC_HAS_RESTRICT)

# clock_gettime() moved from librt to libc in glibc 2.17 (2012-12-25),
# but for backwards compat, still look inside librt.
check_library_exists(rt clock_gettime "time.h" clock_gettime_in_librt)
if(clock_gettime_in_librt)
    set(CMAKE_REQUIRED_LIBRARIES "-lrt")
endif()
check_c_source_runs("
#include <unistd.h>
/* _POSIX_MONOTONIC_CLOCK is one of: -1, 0, or 200809.
 * -1 indicates the facility isn't reliable and shouldn't be trusted. */
#if !defined(_POSIX_MONOTONIC_CLOCK) || _POSIX_MONOTONIC_CLOCK < 0
#error _POSIX_MONOTONIC_CLOCK missing/invalid
#endif

#include <time.h>
int main(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return 0;
}" JEMALLOC_HAVE_CLOCK_MONOTONIC)

check_c_source_runs("
#include <time.h>
int main(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC_COARSE, &ts);
    return 0;
}" JEMALLOC_HAVE_CLOCK_MONOTONIC_COARSE)
set(CMAKE_REQUIRED_FLAGS)
set(CMAKE_REQUIRED_LIBRARIES)

check_c_source_compiles("
#include <mach/mach_time.h>
int main(void) {
    mach_absolute_time();
    return 0;
}" JEMALLOC_HAVE_MACH_ABSOLUTE_TIME)

check_c_source_runs("
#include <assert.h>
#include <strings.h>
int main() {
    int rv = __builtin_ffsl(0x08);
    assert(rv == 4);
    return 0;
}" HAS_BUILTIN_FFSL)

if(HAS_BUILTIN_FFSL)
    set(JEMALLOC_INTERNAL_FFSLL __builtin_ffsll)
    set(JEMALLOC_INTERNAL_FFSL __builtin_ffsl)
    set(JEMALLOC_INTERNAL_FFS __builtin_ffs)
else()
    check_c_source_runs("
    #include <assert.h>
    #include <strings.h>
    int main() {
        int rv = ffsl(0x08);
        assert(rv == 4);
        return 0);
    }" HAS_FFSL)

    if(HAS_FFSL)
        set(JEMALLOC_INTERNAL_FFSLL "ffsll")
        set(JEMALLOC_INTERNAL_FFSL "ffsl")
        set(JEMALLOC_INTERNAL_FFS "ffs")
    else()
        message(FATAL_ERROR "Neither __builtin_ffs() or ffs() found!")
    endif()
endif()
set(CMAKE_REQUIRED_FLAGS)

check_c_source_compiles("
int main(void) {
    __builtin_popcount(0x08);
    return 0;
}" JEMALLOC_INTERNAL_POPCOUNT)
if(JEMALLOC_INTERNAL_POPCOUNT)
    set(JEMALLOC_INTERNAL_POPCOUNT __builtin_popcount)
endif()

check_c_source_compiles("
int main(void) {
    __builtin_popcountl(0x08);
    return 0;
}" JEMALLOC_INTERNAL_POPCOUNTL)
if(JEMALLOC_INTERNAL_POPCOUNTL)
    set(JEMALLOC_INTERNAL_POPCOUNTL __builtin_popcountl)
endif()

set(CMAKE_REQUIRED_FLAGS "-D_GNU_SOURCE")
check_symbol_exists(memalign "stdlib.h" JEMALLOC_OVERRIDE_MEMALIGN)
check_symbol_exists(valloc "stdlib.h" JEMALLOC_OVERRIDE_VALLOC)
check_symbol_exists(secure_getenv "stdlib.h" JEMALLOC_HAVE_SECURE_GETENV)
check_symbol_exists(sched_getcpu "sched.h" JEMALLOC_HAVE_SCHED_GETCPU)
check_symbol_exists(sched_setaffinity "sched.h" JEMALLOC_HAVE_SCHED_SETAFFINITY)
check_symbol_exists(issetugid "unistd.h" JEMALLOC_HAVE_ISSETUGID)
set(CMAKE_REQUIRED_FLAGS)

check_function_exists(_malloc_thread_cleanup JEMALLOC_MALLOC_THREAD_CLEANUP)
if(JEMALLOC_MALLOC_THREAD_CLEANUP)
    set(force_tls YES)
    list(APPEND wrap_syms _malloc_thread_cleanup)
endif()

check_function_exists(_pthread_mutex_init_calloc_cb JEMALLOC_MUTEX_INIT_CB)
if(JEMALLOC_MUTEX_INIT_CB)
    list(APPEND wrap_syms _malloc_prefork _malloc_postfork)
endif()

if(NOT SBRK_DEPRECATED)
    set(JEMALLOC_DSS 1)
endif()

# reset arguments we augmented previously
set(CMAKE_REQUIRED_FLAGS)

if(JEMALLOC_UTRACE)
    check_c_source_compiles("
#include <sys/types.h>
#include <sys/param.h>
#include <sys/time.h>
#include <sys/uio.h>
#include <sys/ktrace.h>

int main(void) {
    utrace((void *)0, 0);
    return 0;
}" JEMALLOC_UTRACE)
endif()

function(checkMadviseFeature feature var)
check_c_source_runs("
#include <sys/mman.h>
int main(void) {
    madvise((void *)0, 0, ${feature});
    return 0;
}" HAS_MADVISE_${feature})
    if(HAS_MADVISE_${feature})
        set(${var} 1 PARENT_SCOPE)
    endif()
endfunction()

checkMadviseFeature(0 JEMALLOC_HAVE_MADVISE)
checkMadviseFeature(MADV_FREE JEMALLOC_PURGE_MADVISE_FREE)
if(NOT JEMALLOC_PURGE_MADVISE_FREE AND JEMALLOC_HAVE_MADVISE)
    if((CMAKE_SYSTEM_PROCESSOR MATCHES "i686" OR
        CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64") AND
        CMAKE_SYSTEM_NAME MATCHES "Linux")
        set(JEMALLOC_PURGE_MADVISE_FREE 1)
        set(JEMALLOC_DEFINE_MADVISE_FREE 1)
    endif()
endif()

checkMadviseFeature(MADV_DONTNEED JEMALLOC_PURGE_MADVISE_DONTNEED)
checkMadviseFeature(MADV_DONTDUMP _dontdump)
checkMadviseFeature(MADV_DODUMP _dodump)
if(_dontdump AND _dodump)
    set(JEMALLOC_MADVISE_DONTDUMP 1)
endif()

checkMadviseFeature(MADV_HUGEPAGE _hugepage)
checkMadviseFeature(MADV_NOHUGEPAGE _nohugepage)
if(_hugepage AND _nohugepage)
    set(JEMALLOC_HAVE_MADVISE_HUGE 1)
endif()

check_c_source_compiles("
int main(void) {
{
    unsigned x = 0;
    int y = __builtin_clz(x);
}
{
    unsigned long x = 0;
    int y = __builtin_clzl(x);
}
    return 0;
}" JEMALLOC_HAVE_BUILTIN_CLZ)

check_c_source_compiles("
#include <os/lock.h>
#include <AvailabilityMacros.h>
int main(void) {
    #if MAC_OS_X_VERSION_MIN_REQUIRED < 101200
    #error not supported
    #else
    os_unfair_lock lock = OS_UNFAIR_LOCK_INIT;
    os_unfair_lock_lock(&lock);
    os_unfair_lock_unlock(&lock);
    #endif
    return 0;
}" JEMALLOC_OS_UNFAIR_LOCK)

check_c_source_compiles("
#include <stddef.h>

extern void (* __free_hook)(void *ptr);
extern void *(* __malloc_hook)(size_t size);
extern void *(* __realloc_hook)(void *ptr, size_t size);

int main(void) {
    void *ptr = 0L;
    if (__malloc_hook) ptr = __malloc_hook(1);
    if (__realloc_hook) ptr = __realloc_hook(ptr, 2);
    if (__free_hook && ptr) __free_hook(ptr);
    return 0;
}" JEMALLOC_GLIBC_MALLOC_HOOK)
# The wrap_syms append for the malloc hooks are in ../CMakeLists.txt because
# it depends on whether or not JEMALLOC_PREFIX is enabled, which we don't
# know at this point in the configuration.
# Potential refactor: move user configuration parsing before compiler detection.

check_c_source_compiles("
#include <stddef.h>

extern void *(* __memalign_hook)(size_t alignment, size_t size);
int main(void) {
   void *ptr = 0L;
    if (__memalign_hook) {
        ptr = __memalign_hook(16, 7);
    }
    return 0;
}" has_memalign_hook)

set(CMAKE_REQUIRED_LIBRARIES "-lpthread")
check_c_source_compiles("
#include <pthread.h>

int main(void) {
    pthread_atfork((void *)0, (void *)0, (void *)0);
    return 0;
}" JEMALLOC_HAVE_PTHREAD_ATFORK)

set(CMAKE_REQUIRED_FLAGS "-D_GNU_SOURCE")
# pthread_setname_np() on macOS only has one argument (the name),
# so even though the function exists on macOS, it fails due to having
# an unexpected prototype.
check_c_source_compiles("
#include <pthread.h>

int main(void) {
    pthread_setname_np(pthread_self(), \"setname_test\");
    return 0;
}" JEMALLOC_HAVE_PTHREAD_SETNAME_NP)
set(CMAKE_REQUIRED_FLAGS)

check_c_source_compiles("
#include <pthread.h>

int main(void) {
    pthread_mutexattr_t attr;
    pthread_mutexattr_init(&attr);
    pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_ADAPTIVE_NP);
    pthread_mutexattr_destroy(&attr);
    return 0;
}" JEMALLOC_HAVE_PTHREAD_MUTEX_ADAPTIVE_NP)

set(CMAKE_REQUIRED_FLAGS "-D_GNU_SOURCE -Werror")
check_c_source_compiles("
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
    char *buffer = (char *)malloc(100);
    char *error = strerror_r(EINVAL, buffer, 100);
    return 0;
}" JEMALLOC_STRERROR_R_RETURNS_CHAR_WITH_GNU_SOURCE)
set(CMAKE_REQUIRED_FLAGS)
set(CMAKE_REQUIRED_LIBRARIES)

check_c_source_compiles("
#include <libkern/OSAtomic.h>
#include <inttypes.h>

int main(void) {
    {
        int32_t x32 = 0;
        volatile int32_t *x32p = &x32;
        OSAtomicAdd32(1, x32p);
    }
    {
        int64_t x64 = 0;
        volatile int64_t *x64p = &x64;
        OSAtomicAdd64(1, x64p);
    }

    return 0;
}" JEMALLOC_OSATOMIC)

check_c_source_runs("
#include <assert.h>
int main(void) {
    int x = 0;
    int val = 1;
    int y = __atomic_fetch_add(&x, val, __ATOMIC_RELAXED);
    int after_add = x;
    assert(after_add == 1);
    return 0;
}" JEMALLOC_GCC_ATOMIC_ATOMICS)

if(force_tls)
    check_c_source_compiles("
    __thread int x;
    int main(void) {
        x = 42;
        return 0;
    }" JEMALLOC_TLS)
endif()

check_c_source_runs("
#include <assert.h>
#include <stdint.h>
int main(void) {
    uint8_t x = 0;
    int32_t val = 1;
    int32_t y = __atomic_fetch_add(&x, val, __ATOMIC_RELAXED);
    int32_t after_add = (int)x;
    assert(after_add == 1);
    return 0;
}" JEMALLOC_GCC_U8_ATOMIC_ATOMICS)

check_c_source_runs("
#include <assert.h>
int main(void) {
    int x = 0;
    int before_add = __sync_fetch_and_add(&x, 1);
    int after_add = x;
    assert(before_add == 0);
    assert(after_add == 1);
    return 0;
}" JEMALLOC_GCC_SYNC_ATOMICS)

check_c_source_runs("
#include <assert.h>
#include <stdint.h>
int main(void) {
    uint8_t x = 0;
    int32_t before_add = __sync_fetch_and_add(&x, 1);
    int32_t after_add = (int)x;
    assert(before_add == 0);
    assert(after_add == 1);
    return 0;
}" JEMALLOC_GCC_U8_SYNC_ATOMICS)

check_c_source_compiles("
void foo(void) {
  __builtin_unreachable();
}

int main(void) {
    return 0;
}" has__builtin_unreachable)
if(has__builtin_unreachable)
    set(JEMALLOC_INTERNAL_UNREACHABLE __builtin_unreachable)
else()
    set(JEMALLOC_INTERNAL_UNREACHABLE abort)
endif()

set(CMAKE_REQUIRED_FLAGS "-Werror -Wall -Wextra -std=c11")
check_c_source_runs("
#include <stdint.h>
#include <assert.h>
#if(__STDC_VERSION__ >= 201112L) && !defined(__STDC_NO_ATOMICS__)
#include <stdatomic.h>
#else
#error Atomics not available
#endif
int main() {
    uint64_t q = 0;
    uint64_t *p = &q;
    uint64_t x = 1;
    volatile atomic_uint_least64_t *a = (volatile atomic_uint_least64_t *)p;
    uint64_t r = atomic_fetch_add(a, x) + x;
    assert(r == 1);
    return 0;
}
" JEMALLOC_C11_ATOMICS)
set(CMAKE_REQUIRED_FLAGS)

UtilCheckTypeSizeValid("void *" LG_SIZEOF_PTR 8 4)
UtilCheckTypeSizeValid("int" LG_SIZEOF_INT 8 4)
UtilCheckTypeSizeValid("long" LG_SIZEOF_LONG 8 4)
UtilCheckTypeSizeValid("long long" LG_SIZEOF_LONG_LONG 8 4)
UtilCheckTypeSizeValid("intmax_t" LG_SIZEOF_INTMAX_T 16 8 4)

# Enable background threads if possible
if(JEMALLOC_HAVE_PTHREAD AND NOT JEMALLOC_OS_UNFAIR_LOCK)
    set(JEMALLOC_BACKGROUND_THREAD 1)
endif()


