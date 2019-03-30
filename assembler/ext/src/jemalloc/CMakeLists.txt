cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

project(jemalloc C CXX)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# Order of includes matters here. Modify at your own risk.
include(${PROJECT_SOURCE_DIR}/build-aux/Utilities.cmake)
include(${PROJECT_SOURCE_DIR}/build-aux/GetSystemSizes.cmake)
include(${PROJECT_SOURCE_DIR}/build-aux/GetSystemAddrBits.cmake)
include(${PROJECT_SOURCE_DIR}/build-aux/DetectOSFeatures.cmake) # 1
include(${PROJECT_SOURCE_DIR}/build-aux/DetectCompilerFeatures.cmake) # 2

# Tell CMake where to find default headers:
include_directories(${PROJECT_SOURCE_DIR}/include)

# Also tell CMake where to find out-of-source generated headers:
include_directories(${PROJECT_BINARY_DIR}/include)
include_directories(${PROJECT_BINARY_DIR}/include/jemalloc)

# User configuration options are maintained in an include for cleanliness:
include(${PROJECT_SOURCE_DIR}/build-aux/UserCompileOptions.cmake)

# Enable smallocx by default too
set(JEMALLOC_EXPERIMENTAL_SMALLOCX_API 1)

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

# If user's compiler supports compileFlag, add it to our build settings
function(setCompileFlag compileFlag)
    # We have to remove all dashes and equals because the test
    # adds defines of -D${flag}, so obviously when a flag has dashes and
    # equals the set won't work (like #define -std=gnu11).
    # So we use the built-in C-conforming generator to fix our input.
    #
    # Sadly, this makes some configure output ugly because CMake doesn't
    # translate anything until build time (also check_c_compiler_flag doesn't
    # have an option for an alternative pretty print string).
    # Alternatively, manually replace bad characters with with:
    #    string(REGEX REPLACE "^-|=|-" "" removeDashEQ ${compileFlag})
    string(MAKE_C_IDENTIFIER ${compileFlag} removeDashEQ)
    check_c_compiler_flag(${compileFlag} ${removeDashEQ})

    if(${removeDashEQ})
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${compileFlag}" PARENT_SCOPE)
    endif()
endfunction()

setCompileFlag("-std=gnu11")

setCompileFlag("-Wall")
setCompileFlag("-Wextra")
setCompileFlag("-Wshorten-64-to-32")
setCompileFlag("-Wsign-compare")
setCompileFlag("-Wundef")
setCompileFlag("-Wno-format-zero-length")

setCompileFlag("-pipe")
setCompileFlag("-g3")

setCompileFlag("-O3")
setCompileFlag("-funroll-loops")

add_definitions(-D_GNU_SOURCE)

# If we can't be gnu11 for some reason, try gnu99 instead
if(NOT (${CMAKE_C_FLAGS} MATCHES "std=gnu11"))
    setCompileFlag("-std=gnu99")
endif()

# If we are gnu11 or gnu99, we know we can use restrict/__restrict__
if((${CMAKE_C_FLAGS} MATCHES "std=gnu11") OR
   (${CMAKE_C_FLAGS} MATCHES "std=gnu99"))
    set(JEMALLOC_HAS_RESTRICT 1)
endif()

# Ask clang to not complain about clearly okay things
if(CMAKE_C_COMPILER_ID MATCHES "Clang")
    setCompileFlag("-Wno-tautological-pointer-compare")
endif()

# Check for C++14 support too
# Original comment from configure.ac about C++ support:
#   Require at least c++14, which is the first version to support sized
#   deallocation.  C++ support is not compiled otherwise.
check_cxx_compiler_flag("-std=gnu++14" HAS_CXX14)
if(HAS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14")
else()
    check_cxx_compiler_flag("-std=c++14" HAS_CXX14)
    if(HAS_CXX14)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
    else()
        message(STATUS "C++14 compiler not found. Not building C++ wrapper.")

        # C++14 support not available. Disable C++ support.
        set(JEMALLOC_CXX_DISABLE TRUE CACHE INTERNAL "Disable CXX")
    endif()
endif()

if(JEMALLOC_VERSION_USER)
    # If user requested specific version number, use directly
    set(jemalloc_version ${JEMALLOC_VERSION_USER})
else()
    # else, introspect the current git repo to find numbers and hashes
    set(jemalloc_version "0.0.0-0-g0000000000000000000000000000000000000000")
    set(jemalloc_version_major 0)
    set(jemalloc_version_minor 0)
    set(jemalloc_version_bugfix 0)
    set(jemalloc_version_nrev 0)
    set(jemalloc_version_gid "0")

    find_package(Git)

    # GetAndParseVersion() populates jemalloc_version when it returns
    GetAndParseVersion()
endif()

if(with-mangling)
    # We are expecting entries separated by a comma
    # with individual entries split by a ':' as in n:m
    # Convert that into a CMake list of ';' separated pairs
    string(REPLACE "," ";" MANGLING_MAP  ${with-mangling})
endif()

# JEMALLOC_PREFIX needs two forms for building:
# JEMALLOC_PREFIX itself is a quoted string: "\"je_\""
# but ALSO the headers use 'je_' as a non-quoted string for
# appending to other symbols directly (and we also use it for
# this purpose and for generating prefixed filenames)
if(JEMALLOC_PREFIX)
    string(REPLACE "\"" "" je_ ${JEMALLOC_PREFIX})
else()
    set(je_ "je_")
endif()

include(CheckFunctionExists)
# function: a function name
# collection: if `function` exists, list for appending function names
# enableOverride: variable created when `function` exists
macro(checkFunctionAppendSyms function collection enableOverride)
    check_function_exists(${function} exists-${function})
    if(exists-${function})
        list(APPEND ${collection} ${function})
        set(${enableOverride} 1)
    endif()
endmacro()

if(NOT JEMALLOC_PREFIX)
    # if no prefix, wrap native allocator functions
    checkFunctionAppendSyms(__libc_calloc wrap_syms
        JEMALLOC_OVERRIDE___LIBC_CALLOC)
    checkFunctionAppendSyms(__libc_free wrap_syms
        JEMALLOC_OVERRIDE___LIBC_FREE)
    checkFunctionAppendSyms(__libc_malloc wrap_syms
        JEMALLOC_OVERRIDE___LIBC_MALLOC)
    checkFunctionAppendSyms(__libc_memalign wrap_syms
        JEMALLOC_OVERRIDE___LIBC_MEMALIGN)
    checkFunctionAppendSyms(__libc_realloc wrap_syms
        JEMALLOC_OVERRIDE___LIBC_REALLOC)
    checkFunctionAppendSyms(__libc_valloc wrap_syms
        JEMALLOC_OVERRIDE___LIBC_VALLOC)
    checkFunctionAppendSyms(__libc_posix_memalign wrap_syms
        JEMALLOC_OVERRIDE___POSIX_MEMALIGN)

    if(JEMALLOC_GLIBC_MALLOC_HOOK)
        list(APPEND wrap_syms __free_hook __malloc_hook __realloc_hook)
    endif()

    if(has_memalign_hook AND NOT JEMALLOC_PREFIX)
        set(JEMALLOC_GLIBC_MEMALIGN_HOOK 1)
        list(APPEND wrap_syms __memalign_hook)
    endif()
endif()

set(JEMALLOC_PREFIX_JET "jet_")

# Base 2 log of minimum tiny size class to support
set(LG_TINY_MIN 3)
if(with-lg-tiny-min)
    set(LG_TINY_MIN ${with-lg-tiny-min})
endif()

# Base 2 log of minimum allocation alignment
set(LG_QUANTA "3 4")

#if(with-lg-quantum)
# Convert to a CMake list
#string(REPLACE "," ";" LG_QUANTA  ${with-lg-quantum})
#  set(LG_QUANTA ${with-lg-quantum})
#set(LG_QUANTUM ${LG_QUANTA})
#endif()

##############################################################################
# On linux, can detect by parsing :
# hpsk=`cat /proc/meminfo 2>/dev/null | \
#          grep -e '^Hugepagesize:[[:space:]]\+[0-9]\+[[:space:]]kB$' | \
#          awk '{print $2}'`
# then adding 1 to the result
set(LG_HUGEPAGE 21)

set(PUBLIC_SYM
    aligned_alloc
    calloc
    dallocx
    free
    mallctl
    mallctlbymib
    mallctlnametomib
    malloc
    malloc_conf
    malloc_message
    malloc_stats_print
    malloc_usable_size
    mallocx
    nallocx
    posix_memalign
    rallocx
    realloc
    sallocx
    sdallocx
    xallocx
    valloc
)

checkFunctionAppendSyms(memalign PUBLIC_SYM JEMALLOC_OVERRIDE_MEMALIGN)
checkFunctionAppendSyms(valloc PUBLIC_SYM JEMALLOC_OVERRIDE_VALLOC)

if(JEMALLOC_SMALLOCX_ENABLE)
    list(APPEND PUBLIC_SYM smallocx_${jemalloc_version_gid})
endif()

# Generate prefixed symbol names for awk nm extractors
foreach(sym ${PUBLIC_SYM})
    list(APPEND PREFIXED_JE  "${je_}${sym}")
    list(APPEND PREFIXED_JET "jet_${sym}")
endforeach()

# system API overrides (these don't get prefixed)
foreach(sym ${wrap_syms})
    list(APPEND PREFIXED_JE  "${sym}")
    list(APPEND PREFIXED_JET "${sym}")
endforeach()

macro(createSymbolsAwk end syms)
add_custom_command(
    OUTPUT
    "${PROJECT_BINARY_DIR}/include/jemalloc/internal/private_symbols${end}.awk"
    COMMAND
    "${PROJECT_SOURCE_DIR}/include/jemalloc/internal/private_symbols.sh"
    "${NM_PREFIX}"
    ${syms}
    >
    "${PROJECT_BINARY_DIR}/include/jemalloc/internal/private_symbols${end}.awk"
    VERBATIM
    # Add a depends on the source file itself because ${syms} is a static list
    # of function names defined here in CMakeLists.txt, so we can't tell the
    # input to this command is newer than the output file or not.
    DEPENDS
    "${PROJECT_SOURCE_DIR}/CMakeLists.txt"
)
endmacro()

createSymbolsAwk("" "${PREFIXED_JE}")
createSymbolsAwk("_jet" "${PREFIXED_JET}")

add_custom_target(jemalloc-private-symbols-awk DEPENDS
    "${PROJECT_BINARY_DIR}/include/jemalloc/internal/private_symbols.awk")

add_custom_target(jemalloc-private-symbols-awk-jet DEPENDS
    "${PROJECT_BINARY_DIR}/include/jemalloc/internal/private_symbols_jet.awk")

set(PUBLIC_SYM_FILE
    "${PROJECT_BINARY_DIR}/include/jemalloc/internal/public_symbols.txt")

GeneratePublicSymbolsList(
    "${PUBLIC_SYM}"
    "${MANGLING_MAP}"
    ${je_}
    "${PUBLIC_SYM_FILE}")

foreach(public_in jemalloc_macros.h
                  jemalloc_defs.h
                  jemalloc_protos.h
                  jemalloc_typedefs.h)
    ConfigureFile("${PROJECT_SOURCE_DIR}/include/jemalloc/${public_in}.in"
                  "${PROJECT_BINARY_DIR}/include/jemalloc/${public_in}" True)
endforeach()

# Generate main public header
set(JEMALLOC_HDR
    "${PROJECT_BINARY_DIR}/include/jemalloc/jemalloc${install_suffix}.h"
    CACHE INTERNAL "jemalloc public header")

# These are headers compiled by jemalloc.sh
set(JEMALLOC_HDR_LIST
    jemalloc_defs.h
    jemalloc_rename.h
    jemalloc_macros.h
    jemalloc_protos.h
    jemalloc_typedefs.h
    jemalloc_mangle.h
)

CreateJemallocHeader("${PUBLIC_SYM_FILE}"
                     "${JEMALLOC_HDR_LIST}"
                     "${JEMALLOC_HDR}")

##############################################################################
## Internal headers generation
set(PUBLIC_NAMESPACE_FILE
    "${PROJECT_BINARY_DIR}/include/jemalloc/internal/public_namespace.h")
PublicNamespace(${PUBLIC_SYM_FILE} "${PUBLIC_NAMESPACE_FILE}")

set(PUBLIC_UNNAMESPACE_FILE
    "${PROJECT_BINARY_DIR}/include/jemalloc/internal/public_unnamespace.h")
PublicUnnamespace(${PUBLIC_SYM_FILE} "${PUBLIC_UNNAMESPACE_FILE}")

# Configure internal headers
ConfigureFile(
    "${PROJECT_SOURCE_DIR}/include/jemalloc/internal/jemalloc_preamble.h.in"
    "${PROJECT_BINARY_DIR}/include/jemalloc/internal/jemalloc_preamble.h" True)

foreach(internal_in jemalloc_internal_defs.h)
    ConfigureFile(
        "${PROJECT_SOURCE_DIR}/include/jemalloc/internal/${internal_in}.in"
        "${PROJECT_BINARY_DIR}/include/jemalloc/internal/${internal_in}" True)
endforeach()

# Headers for test running
ConfigureFile(
    "${PROJECT_SOURCE_DIR}/test/include/test/jemalloc_test_defs.h.in"
    "${PROJECT_BINARY_DIR}/test/include/test/jemalloc_test_defs.h" True)

ConfigureFile(
    "${PROJECT_SOURCE_DIR}/test/include/test/jemalloc_test.h.in"
    "${PROJECT_BINARY_DIR}/test/include/test/jemalloc_test.h" False)

# Use jet_ generate protos_jet
set(je_ "jet_")

# replace prefix only
ConfigureFile(
    "${PROJECT_SOURCE_DIR}/include/jemalloc/jemalloc_protos.h.in"
    "${PROJECT_BINARY_DIR}/include/jemalloc/jemalloc_protos_jet.h" True)

# revert
if(JEMALLOC_PREFIX)
    string(REPLACE "\"" "" je_ ${JEMALLOC_PREFIX})
else()
    set(je_ "je_")
endif()

# Enable testing when in Debug mode
if(CMAKE_BUILD_TYPE MATCHES Debug)
    include_directories(${PROJECT_BINARY_DIR}/test/include)
    include_directories(${PROJECT_BINARY_DIR}/include/jemalloc/internal)

    include(CTest)
    enable_testing()

    # check target is scoped with jemalloc- because when jemalloc
    # is included as cmake dependency in other projects, we don't
    # want to conflict with any other 'check' targets
    add_custom_target(${PROJECT_NAME}-check COMMAND ${CMAKE_CTEST_COMMAND})

    add_subdirectory(test)
endif()

add_subdirectory(src)
