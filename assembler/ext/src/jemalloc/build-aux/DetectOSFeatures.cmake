# This file configures OS and/or hardware support options

include(CheckSymbolExists)

set(JEMALLOC_MAPS_COALESCE 1)

# Adapted from configure.ac section: case "${host}" in
string(TOLOWER ${CMAKE_SYSTEM_NAME} platform) # case insensitive compare
set(NM_PREFIX "")
if(platform MATCHES "darwin")
    set(abi "macho")
    set(SBRK_DEPRECATED YES)
    set(force_tls 0)
    set(JEMALLOC_ZONE 1)
    set(NM_PREFIX "_") # macho symbols are prefixed with _ but elf isn't
elseif(platform MATCHES "linux")
    set(abi "elf")
    set(JEMALLOC_PURGE_MADVISE_DONTNEED_ZEROS 1)
    set(JEMALLOC_HAS_ALLOCA_H 1)
    set(JEMALLOC_PROC_SYS_VM_OVERCOMMIT_MEMORY 1)
    set(JEMALLOC_THREADED_INIT 1)
    set(JEMALLOC_USE_CXX_THROW 1)
elseif(platform MATCHES "freebsd")
    set(abi "elf")
    set(JEMALLOC_SYSCTL_VM_OVERCOMMIT 1)
    set(force_lazy_lock 1)
elseif(platform MATCHES "dragonfly")
    set(abi "elf")
elseif(platform MATCHES "openbsd")
    set(abi "elf")
    set(force_tls 0)
elseif(platform MATCHES "linux-android")
    set(abi "elf")
    set(JEMALLOC_PURGE_MADVISE_DONTNEED_ZEROS 1)
    set(JEMALLOC_HAS_ALLOCA_H 1)
    set(JEMALLOC_PROC_SYS_VM_OVERCOMMIT_MEMORY 1)
    set(JEMALLOC_THREADED_INIT 1)
    set(JEMALLOC_C11_ATOMICS 1)
elseif(platform MATCHES "kfreebsd")
    set(abi "elf")
    set(JEMALLOC_HAS_ALLOCA_H 1)
    set(JEMALLOC_SYSCTL_VM_OVERCOMMIT 1)
    set(JEMALLOC_THREADED_INIT 1)
    set(JEMALLOC_USE_CXX_THROW 1)
elseif(platform MATCHES "netbsd")
    check_symbol_exists("__ELF__" "stdlib.h" useELF)
    if (useELF)
        set(abi "elf")
    else()
        set(abi "aout")
    endif()
elseif(platform MATCHES "solaris2")
    set(abi "elf")
elseif(platform MATCHES "ibm-aix")
    set(abi "xcoff")
elseif(platform MATCHES "mingw" OR platform MATCHES "cygwin")
    set(JEMALLOC_MAPS_COALESCE 0)
    set(abi "pecoff")
else()
    message(WARNING "Unsupported operating system: ${platform}")
    set(abi "elf")
endif()

# Only run check if we haven't run it before
if(CMAKE_SYSTEM_PROCESSOR MATCHES "i686" OR
   CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
    check_c_source_compiles("
int main(void) {
    __asm__ volatile(\"pause\");
    return 0;
    }" HAVE_CPU_SPINWAIT)
    if(HAVE_CPU_SPINWAIT)
        set(CPU_SPINWAIT "__asm__ volatile(\"pause\")" CACHE INTERNAL "CS")
    endif()
elseif(MSVC)
    check_c_source_compiles("
int main(void) {
    _mm_pause();
    return 0;
}" HAVE_CPU_SPINWAIT)
    if(HAVE_CPU_SPINWAIT)
        set(CPU_SPINWAIT "_mm_pause()" CACHE INTERNAL "CS")
    else()
        check_c_source_compiles("
#include <Windows.h>
int main(void) {
    YieldProcessor();
    return 0;
}" HAVE_CPU_SPINWAIT)
        if(HAVE_CPU_SPINWAIT)
            set(CPU_SPINWAIT "YieldProcessor()" CACHE INTERNAL "CS")
        endif()
    endif()
endif()

# If defined, use munmap() to unmap freed chunks, rather than storing them for
# later reuse.  This is disabled by default on Linux because common sequences
# of mmap()/munmap() calls will cause virtual memory map holes.
# But it is enabled by default on Windows
if(NOT UNIX)
    set(JEMALLOC_MUNMAP 1)
endif()

if(NOT LG_VADDR)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64") # maybe "arm64" too?
        set(LG_VADDR 48)
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
        GetSystemAddrBits(LG_VADDR)
        # Cache result so we don't run the check every time
        set(LG_VADDR ${LG_VADDR} CACHE INTERNAL "System Address Bits")
    else()
        message(FATAL_ERROR "Unknown number of VADDR bits")
    endif()
endif()

if((NOT LG_PAGE) OR (NOT LG_CACHELINE))
    if(NOT SYSTEM_PAGE_SIZE)
        GetSystemSizes(SYSTEM_PAGE_SIZE SYSTEM_CACHE_SIZE)
        set(SYSTEM_PAGE_SIZE ${SYSTEM_PAGE_SIZE} CACHE INTERNAL "Page Size")
        set(SYSTEM_CACHE_SIZE ${SYSTEM_CACHE_SIZE} CACHE INTERNAL "Cache Line Size")
    endif()
    lg(${SYSTEM_PAGE_SIZE} LG_PAGE)
    lg(${SYSTEM_CACHE_SIZE} LG_CACHELINE)
    set(LG_PAGE ${LG_PAGE} CACHE INTERNAL "lg page size")
    set(LG_CACHELINE ${LG_CACHELINE} CACHE INTERNAL "lg cache line")
    set(CACHELINE ${SYSTEM_CACHE_SIZE} CACHE INTERNAL "alias for macro replacement")
endif()

if(NOT WIN32)
    include(CheckCSourceCompiles)

    # Check if syscall(2) is usable.
    # Treat warnings as errors, so that e.g. OS X 10.12's warning prevents use.
    set(CMAKE_REQUIRED_FLAGS "-Werror -Wall")
    check_c_source_compiles("
    #include <sys/syscall.h>
    #include <unistd.h>
    int main(void) {
        syscall(SYS_write, 2, \"hello\", 5);
        return 0;
    }" JEMALLOC_USE_SYSCALL)
    set(CMAKE_REQUIRED_FLAGS)
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "FreeBSD")
    # Relevant for FreeBSD only
    set(JEMALLOC_ATOMIC9 1)
endif()

# OS X / iOS / Darwin settings
if(APPLE)
    set(LG_HUGEPAGE 21)
endif()
