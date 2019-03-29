# We need this as a SOURCES array and not _just_ an OBJECT library because
# these sources will be compiled multiple times (with name-prefix jet-* and
# without (regular distribution)).
set(JEMALLOC_CMAKE_SOURCES
    arena.c
    background_thread.c
    base.c
    bin.c
    bitmap.c
    ckh.c
    ctl.c
    div.c
    extent.c
    extent_dss.c
    extent_mmap.c
    hash.c
    hook.c
    jemalloc.c
    large.c
    log.c
    malloc_io.c
    mutex.c
    mutex_pool.c
    nstime.c
    pages.c
    prng.c
    prof.c
    rtree.c
    sc.c
    stats.c
    sz.c
    tcache.c
    test_hooks.c
    ticker.c
    tsd.c
    witness.c
    )

if(APPLE)
    list(APPEND JEMALLOC_CMAKE_SOURCES zone.c)
endif()

# These objects only have code in debug builds (due to macros)
if(CMAKE_BUILD_TYPE MATCHES Debug)
    list(APPEND JEMALLOC_CMAKE_SOURCES
        hash.c
        prng.c
        ticker.c)
    # No symbols even in debug mode
    # mb.c
endif()

if(NOT JEMALLOC_CXX_DISABLE)
    # C++14 support exists, so add C++ file and configure C++ options
    list(APPEND JEMALLOC_CMAKE_SOURCES
        jemalloc_cpp.cpp)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g3 -O3")
endif()

if(enable_valgrind)
    list(APPEND C_SRCS valgrind.c)
endif()

# fix mac rpath complaints from cmake
if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif()

find_program(JEMALLOC_AWK_EXEC awk)
function(ExtractSymbolsGenerateNamespaceHeader target symExt awk dir hdr)
# Don't compile -syms with optimization so it builds faster.
# The target is only for symbol extraction so we don't care about performance.
target_compile_options(${target} PRIVATE "-O0")
target_compile_definitions(${target} PRIVATE
    JEMALLOC_NO_PRIVATE_NAMESPACE
    _GNU_SOURCE)

# Stash the nm symbol extraction output in a subdirectory so we don't
# clutter build/src itself with dozens of temporary files.
file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/src/syms/")

# Generate per-file output: nm -a $obj | awk -f private_symbols.awk > $obj.sym
foreach(sym ${JEMALLOC_CMAKE_SOURCES})
    get_filename_component(filename ${sym} NAME)
    set(obj ${filename}.o)
    set(objectDir "${PROJECT_BINARY_DIR}/src/CMakeFiles/${dir}")
    add_custom_command(
        TARGET
        ${target}
        PRE_LINK
        BYPRODUCTS
        "${PROJECT_BINARY_DIR}/src/syms/${obj}.${symExt}"
        COMMAND
        ${CMAKE_NM} -a "${objectDir}/${obj}" |
        ${JEMALLOC_AWK_EXEC} -f
        "${PROJECT_BINARY_DIR}/include/jemalloc/internal/${awk}"
        >
        "${PROJECT_BINARY_DIR}/src/syms/${obj}.${symExt}"
        COMMENT
        "Extracting symbols from ${obj}"
)
endforeach()

add_custom_command(
    TARGET
    ${target}
    PRE_LINK
    BYPRODUCTS
    "${PROJECT_BINARY_DIR}/include/jemalloc/internal/${hdr}"
    COMMAND
    "${PROJECT_SOURCE_DIR}/include/jemalloc/internal/private_namespace.sh"
    "${PROJECT_BINARY_DIR}/src/syms/*.${symExt}"
    >
    "${PROJECT_BINARY_DIR}/include/jemalloc/internal/${hdr}"
    COMMENT
    "Generating ${target} ${hdr}"
)

# This is a hack because our generated headers are created _after_
# the generated objects, causing the objects to be OLDER than the
# included header they depend on, which causes a never-ending
# rebuild loop of all sources.
# The simple hack here is to just go back and touch all the object
# files after the common header is generated but before the final
# archives are linked.
foreach(sym ${JEMALLOC_CMAKE_SOURCES})
    get_filename_component(filename ${sym} NAME)
    set(obj ${filename}.o)
    set(objectDir "${PROJECT_BINARY_DIR}/src/CMakeFiles/${dir}")
    add_custom_command(
        TARGET
        ${target}
        PRE_LINK
        COMMAND
        ${CMAKE_COMMAND} -E touch "${objectDir}/${obj}"
        #        COMMENT
        #        "Touching ${target} ${objectDir}/${obj}"
    )
endforeach()

endfunction()

add_library(${PROJECT_NAME}-syms STATIC ${JEMALLOC_CMAKE_SOURCES})

# Stash "syms" archive in a subdirectory so it doesn't confuse users
set_target_properties(${PROJECT_NAME}-syms PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "syms/")

# Tell CMake we need the awk template generated before running symbol extraction
add_dependencies(${PROJECT_NAME}-syms jemalloc-private-symbols-awk)

ExtractSymbolsGenerateNamespaceHeader(
    ${PROJECT_NAME}-syms
    "sym"
    "private_symbols.awk"
    "${PROJECT_NAME}-syms.dir"
    "private_namespace.h"
)

# Bring sources into an object library so we can include elsewhere.
# (also so we don't compile twice for both -static and -shared)
add_library(${PROJECT_NAME} OBJECT ${JEMALLOC_CMAKE_SOURCES})

# The original library, delivery product
add_library(${PROJECT_NAME}-static STATIC $<TARGET_OBJECTS:${PROJECT_NAME}>)

# Now add shared library. Needed for integration tests and benchmarks
add_library(${PROJECT_NAME}-shared SHARED $<TARGET_OBJECTS:${PROJECT_NAME}>)

add_dependencies(${PROJECT_NAME} ${PROJECT_NAME}-syms)

# Set the name of static and shared lib filenames to be "libjemalloc.*"
# instead of their full CMake libary names of libjemalloc-static.a, etc.
set_target_properties(${PROJECT_NAME}-shared PROPERTIES
    OUTPUT_NAME ${PROJECT_NAME})
set_target_properties(${PROJECT_NAME}-static PROPERTIES
    OUTPUT_NAME ${PROJECT_NAME})

# Need to export from DLL
target_compile_definitions(${PROJECT_NAME}-shared PRIVATE DLLEXPORT)

if(with-jemalloc-prefix)
    target_compile_definitions(${PROJECT_NAME}-static PRIVATE JEMALLOC_MANGLE)
    target_compile_definitions(${PROJECT_NAME}-shared PRIVATE JEMALLOC_MANGLE)
endif()

if(enable-prefixed-jemalloc-library)
    add_library(${PROJECT_NAME}-jet OBJECT ${JEMALLOC_CMAKE_SOURCES})
    add_library(${PROJECT_NAME}-jet-static STATIC
        $<TARGET_OBJECTS:${PROJECT_NAME}-jet>)

    target_compile_definitions(${PROJECT_NAME}-jet
        PRIVATE
        "${C_JETLIB_DEFS}"
    )

    add_library(${PROJECT_NAME}-jet-syms
        STATIC
        ${JEMALLOC_CMAKE_SOURCES}
    )

    target_compile_definitions(${PROJECT_NAME}-jet-syms
        PRIVATE
        "${C_JETLIB_DEFS}"
    )

    set_target_properties(${PROJECT_NAME}-jet-syms PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "syms/")

    add_dependencies(${PROJECT_NAME}-jet-syms
        jemalloc-private-symbols-awk-jet)

    add_dependencies(${PROJECT_NAME}-jet-static
        ${PROJECT_NAME}-jet-syms)

    add_dependencies(${PROJECT_NAME}-jet ${PROJECT_NAME}-jet-syms)

    ExtractSymbolsGenerateNamespaceHeader(
        ${PROJECT_NAME}-jet-syms
        "sym.jet"
        "private_symbols_jet.awk"
        "${PROJECT_NAME}-jet-syms.dir"
        "private_namespace_jet.h"
    )
endif()
