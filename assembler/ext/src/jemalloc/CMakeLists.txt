cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
# Require CMake >= 3.0 for check_c_compiler_flag support

project(jemalloc C CXX)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CTest)

include(${PROJECT_SOURCE_DIR}/build-aux/Utilities.cmake)
include(${PROJECT_SOURCE_DIR}/build-aux/GetSystemSizes.cmake)
include(${PROJECT_SOURCE_DIR}/build-aux/GetSystemAddrBits.cmake)
include(${PROJECT_SOURCE_DIR}/build-aux/DetectOSFeatures.cmake) # 1
include(${PROJECT_SOURCE_DIR}/build-aux/DetectCompilerFeatures.cmake) # 2

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_BINARY_DIR}/include)
include_directories(${PROJECT_BINARY_DIR}/include/jemalloc)

# User configuration options are maintained in an include for cleanliness:
include(${PROJECT_SOURCE_DIR}/build-aux/UserCompileOptions.cmake)

# Enable smallocx by default too
set(JEMALLOC_EXPERIMENTAL_SMALLOCX_API 1)

# install_prefix - installation directory prefix
# with-xslroot=<path>  XSL stylesheet root path

# Installation
# Munge install path variables.
# All of these can simply we be passed as
# Cmake arguments
if(NOT install_prefix)
    set(install_prefix "/usr/local")
endif()

# We do not support exec_prefix until we find it is necessary
if(NOT exec_prefix)
    set(exec_prefix ${install_prefix})
endif()

set(PREFIX ${install_prefix})

if(enable-prefixed-jemalloc-library)
    set(with-jemalloc-prefix ON)
endif()

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

# If user's compiler supports compileFlag, add it to our build settings
macro(setCompileFlag compileFlag)
	# We have to remove all dashes and equals because the test
    # addes defines of -D${flag}, so obviously when a flag has dashes and
    # equals the set won't work (like #define -std=gnu11).
    # So we use the built-in C-conforming generator to fix our input.
    # 
    # Sadly, this makes some configure output ugly because CMake doesn't
    # translate anything until build time (also check_c_compiler_flag doesn't
    # have an option for an alternative pretty print string).
    # Alternatively, manually replace bad characters with with:
    #	string(REGEX REPLACE "^-|=|-" "" removeDashEQ ${compileFlag})
    set(removeDashEQ $<MAKE_C_IDENTIFIER:${compileFlag}>)
    check_c_compiler_flag(${compileFlag} ${removeDashEQ})

	if(${removeDashEQ})
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${compileFlag}")
	endif()
endmacro()

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
setCompileFlag("-D_GNU_SOURCE")

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

if(JEMALLOC_VERSION)
    # If user requested specific version number, use directly
    set(jemalloc_version ${JEMALLOC_VERSION})
else()
    # else, introspect the current git repo to find numbers and hashes
    set(jemalloc_version "0.0.0-0-g0000000000000000000000000000000000000000")
    set(jemalloc_version_major 0)
    set(jemalloc_version_minor 0)
    set(jemalloc_version_bugfix 0)
    set(jemalloc_version_nrev 0)
    set(jemalloc_version_gid "0")

    find_package(Git)
    GetAndParseVersion()
endif()

if(with-mangling)
    # We are expecting entries separated by a comma
    # with individual entries split by a ':' as in n:m
    # Convert that into a CMake list of ';' separated pairs
    string(REPLACE "," ";" MANGLING_MAP  ${with-mangling})
endif()

# Set the default API prefix for public
# Need quotes so the preprocessor concats two strings
set(JEMALLOC_PREFIX "\"je_\"")
set(je_ je_) # in-header prefix for je_blahblah

if(with-jemalloc-prefix)
    set(JEMALLOC_PREFIX ${JEMALLOC_PREFIX})
    set(JEMALLOC_PRIVATE_NAMESPACE ${JEMALLOC_PREFIX})
endif()

# Uppercase copy of the JEMALLOC_PREFIX
if(JEMALLOC_PREFIX)
    string(TOUPPER ${JEMALLOC_PREFIX} JEMALLOC_CPREFIX)
endif()

# Disable exporting jemalloc public APIs
# We need to define the var to whitespace string
# as empty strings will not be defined in CMake
# Non-empty definition is necessary so
if(without-export)
    set(JEMALLOC_EXPORT " ")
endif()

if(with-install-suffix)
    set(INSTALL_SUFFIX ${with-install-suffix})
    set(install_suffix ${with-install-suffix})
endif()

# Build debugging code (implies --enable-ivsalloc)
if(enable-debug)
    set(JEMALLOC_DEBUG 1)
endif()

set(JEMALLOC_PREFIX_JET jet_)

# Support xmalloc option
set(JEMALLOC_XMALLOC 0)
if(enable-xmalloc)
    set(JEMALLOC_XMALLOC 1)
endif()

# Support cache-oblivious allocation alignment by default.
# If defined, explicitly attempt to more uniformly distribute large allocation
# pointer alignments across all cache indices.
set(JEMALLOC_CACHE_OBLIVIOUS 1)
if(disable-cache-oblivious)
    set(JEMALLOC_CACHE_OBLIVIOUS 0)
endif()

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

# Base 2 log of system page size

if(with-lg-page)
    set(LG_PAGE ${with-lg-page})
endif()

# Base 2 logs of system page sizes to support
set(LG_PAGE_SIZES ${LG_PAGE})
if(with-lg-page-sizes)
    string(REPLACE "," ";" LG_PAGE_SIZES  ${with-lg-page-sizes})
endif()

# Base 2 log of size classes per doubling
set(LG_SIZE_CLASS_GROUP 2)
if(with-lg-size-class-group)
    set(LG_SIZE_CLASS_GROUP ${with-lg-size-class-group})
endif()

set(JEMALLOC_HAVE_ISSETUGID 0)
set(JEMALLOC_MALLOC_THREAD_CLEANUP 0)

set(JEMALLOC_MUTEX_INIT_CB 0)

##############################################################################
# enable-lazy-lock
set(JEMALLOC_LAZY_LOCK 0)
if(NOT enable-lazy-lock)
    if(force_lazy_lock)
        message(STATUS "Forcing lazy-lock to avoid allocator/threading bootstrap issues")
        set(enable-lazy-lock ON)
    endif()
endif()

if(enable-lazy-lock)
    if(${abi} STREQUAL "pecoff")
        message(STATUS "Forcing no lazy-lock because thread creation monitoring is unimplemented")
        set(enable-lazy-lock OFF)
    endif()
endif()

if(enable-lazy-lock)
endif()

##############################################################################

# On linux, can detect by parsing :
# hpsk=`cat /proc/meminfo 2>/dev/null | \
#          grep -e '^Hugepagesize:[[:space:]]\+[0-9]\+[[:space:]]kB$' | \
#          awk '{print $2}'`
# then adding 1 to the result
set(LG_HUGEPAGE 21)

##############################################################################
# Generate configured public headers for concatenation
# Public Headers in for configuring
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

if(JEMALLOC_SMALLOCX_ENABLE)
    list(APPEND PUBLIC_SYM
        smallocx_${jemalloc_version_gid}
    )
endif()

# Generate prefixed symbol names for awk nm extractors
foreach(sym ${PUBLIC_SYM})
    list(APPEND PREFIXED_JE  "${je_}${sym}")
    list(APPEND PREFIXED_JET "jet_${sym}")
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
    DEPENDS
    "${PROJECT_SOURCE_DIR}/include/jemalloc/internal/private_symbols.sh"
    VERBATIM
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
    "${PROJECT_BINARY_DIR}/include/jemalloc/jemalloc${install_suffix}.h")

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
#ConfigureFile(
#    "${PROJECT_SOURCE_DIR}/test/include/test/jemalloc_test_defs.h.in"
#    "${PROJECT_BINARY_DIR}/test/include/test/jemalloc_test_defs.h" True)

#ConfigureFile(
#    "${PROJECT_SOURCE_DIR}/test/include/test/jemalloc_test.h.in"
#    "${PROJECT_BINARY_DIR}/test/include/test/jemalloc_test.h" False)

# Use jet_ generate protos_jet
set(je_ jet_)

# replace prefix only
ConfigureFile(
    "${PROJECT_SOURCE_DIR}/include/jemalloc/jemalloc_protos.h.in"
    "${PROJECT_BINARY_DIR}/include/jemalloc/jemalloc_protos_jet.h" True)

# revert
set(je_ je_)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(enable-prefixed-jemalloc-library ON)
    include_directories(${PROJECT_BINARY_DIR}/test/include)
endif()

# TODO: all of these variables need to be refactored (here and src/ and test/).
# All these things are from the original jemalloc-cmake repo which had
# some oddities we are still trying to flush out.
if(enable-prefixed-jemalloc-library)
    # jet builds use header in jemalloc/internal binary path
    include_directories(${PROJECT_BINARY_DIR}/include/jemalloc/internal)
    set(C_JETLIB_DEFS "JEMALLOC_JET")

    set(C_UTIL_INTEGRATION_DEFS "")
    set(C_TESTLIB_INTEGRATION_DEFS "JEMALLOC_INTEGRATION_TEST")
    set(C_TESTLIB_INTEGRATION ${PROJECT_NAME}IntegrationTest${install_suffix})

    set(C_TESTLIB_UNIT_DEFS "JEMALLOC_UNIT_TEST")
    set(C_TESTLIB_STRESS_DEF "JEMALLOC_STRESS_TEST;JEMALLOC_STRESS_TESTLIB")
    set(C_TESTLIB_STRESS ${PROJECT_NAME}TestStress${install_suffix})

    set(C_UNITETEST_DEFS "JEMALLOC_UNIT_TEST")

    if(with-jemalloc-prefix)
        set(C_JETLIB_DEFS  "${C_JETLIB_DEFS};JEMALLOC_MANGLE")
        set(C_UTIL_INTEGRATION_DEFS
            "${C_UTIL_INTEGRATION_DEFS};JEMALLOC_MANGLE")
    endif()
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
    enable_testing()

    # check target is scoped with -jemalloc because when jemalloc
    # is included as cmake dependency in other projects, we don't
    # want to conflict with any other 'check' targets
    add_custom_target(check-jemalloc COMMAND ${CMAKE_CTEST_COMMAND})

    add_subdirectory(test)
endif()

add_subdirectory(src)
