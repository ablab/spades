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

    # Also, test binaries link against jemalloc-static,
    # so override any user preference and create static library when testing.
    set(JEMALLOC_BUILD_STATIC ON)

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
function(ExtractSymbolsGenerateNamespaceHeader target top-target symExt awk dir hdr)
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

    # Control symbol rebulding either when:
    # - only rebuild symbols for changed files (makes CMake output cluttered)
    # - rebuild all symbols when any source file changes (nicer CMake output
    #   for end users, but shows each symbol extract line on every source change)
    set(EXTRACT_PER_FILE OFF)

    # There are two ways of generating the symbol files and the
    # above feature toggle EXTRACT_PER_FILE controls which we're attempting:
    # - when ON: create a change-notification target for _each_ source file so
    #   symbols will get updated only for updated source files when source
    #   files change
    # - OR -
    # - when OFF: create a library update command so when any source file changes and
    #   the "jemalloc-syms-internal" library is created again, ALL sources
    #   have their symbols extracted.
    #
    # The per-file-dependency approach benefits from only symbols are extracted
    # for updated files during development, but it causes CMake to generate noisy
    # output because *each symbol output* file also appears as a compile target in
    # CMake output. (There's also a problem with on-demand generating the
    # final private_namespace.h file, where in this per-file-symbol mode, it always
    # runs which makes for also noisy CMake output).
    #
    # The per-library-dependency approach benefits from not having noisy CMake output
    # and also having proper header file generation on-demand, BUT the problem is:
    # during development, any updated source file causes all symbols to be regenerated.
    # The regenerate-all-symbols isn't a problem in practice (it's fast), but it
    # does clutter the build output screen.
    #
    # If anybody can figure out how to make per-file dependencies without cluttering
    # the CMake target output list (and without always re-generating the headers
    # because it can't quite track change dependencies correctly), give it a try.
    #
    # Super alternatively:
    #   - if we remove COMMENT tags from the following commands and targets,
    #     the symbol rebuilds can run whenever they want and users won't be bothered
    #     with noisy output.
    if (EXTRACT_PER_FILE)
        set(symTarget ${target}-${obj}.sym)
        add_custom_command(
            OUTPUT
            "${PROJECT_BINARY_DIR}/src/syms/${obj}.${symExt}"
            COMMAND
            ${CMAKE_NM} -a "${objectDir}/${obj}" |
            ${JEMALLOC_AWK_EXEC} -f
            "${PROJECT_BINARY_DIR}/include/jemalloc/internal/${awk}"
            >
            "${PROJECT_BINARY_DIR}/src/syms/${obj}.${symExt}"
            DEPENDS
            "${objectDir}/${obj}"
            COMMENT
            "Extracting symbols from ${obj} for ${target}"
        )

        # Target can only be created once per project, so regen must be blocked
        if (NOT TARGET ${symTarget})
            # CMake can't use output files directly as targets, so we need to wrap
            # our output file in a custom target, then include the custom target as a
            # dependency of the top level target so everything builds in a correct oder.
            add_custom_target(${symTarget}
                DEPENDS
                "${PROJECT_BINARY_DIR}/src/syms/${obj}.${symExt}"
                ${target}
            )

            add_dependencies(${top-target} ${symTarget})
        endif()
    else()
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
            "Extracting symbols from ${obj} for ${target}"
        )
    endif()
endforeach()

# Note: if using EXTRACT_PER_FILE==ON, you need to change the following to:
#       TARGET ${top-target} POST_BUILD
#       because CMake ends up running header generation before symbol extraction
#       which clearly doesn't work for us.
#       BUT—when you change to POST_BUILD on ${top-target} the header generation
#       will ALWAYS run because ${top-target} is a placeholder and doesn't actually
#       build anything, which confuses CMake, thus the rule always runs.
# When using EXTRACT_PER_FILE==OFF, we use TARGET ${target} and PRE_LINK because
# CMake *will* run header generation after the previous PRE_LINK target above.
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
    "${PROJECT_BINARY_DIR}/include/jemalloc/internal/checkAgainst-${hdr}"
    # Now, do some trickery and *only* create a new header if it is actually
    # different from an existing header. This helps us fix Make's problem of
    # thinking a header is "newer" than all object files and constantly
    # causing full rebuilds even when nothing has changed.
    COMMAND
    ${CMAKE_COMMAND} -E copy_if_different
    "${PROJECT_BINARY_DIR}/include/jemalloc/internal/checkAgainst-${hdr}"
    "${PROJECT_BINARY_DIR}/include/jemalloc/internal/${hdr}"
    # DEPENDS on jemalloc.h because experimental symbol names get assigned
    # on each new git commit, so we need to update
    # these generated headers too after any version updates are detected:
    DEPENDS
    "${PROJECT_BINARY_DIR}/include/jemalloc/jemalloc.h"
    COMMENT
    "Generating ${hdr} for ${target}"
)
endfunction()

# We have two "syms" targets because if we want to have per-file symbol
# rebuilding, we need all .o files to be generated first, so the symbol
# extraction ("-syms") will depend on objects created by -syms-internal.
add_library(${PROJECT_NAME}-syms-internal STATIC ${JEMALLOC_CMAKE_SOURCES})
add_dependencies(${PROJECT_NAME}-syms-internal jemalloc-private-symbols-awk)

# Stash "syms" archive in a subdirectory so it doesn't confuse users
set_target_properties(${PROJECT_NAME}-syms-internal PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "syms/")

add_custom_target(${PROJECT_NAME}-syms DEPENDS ${PROJECT_NAME}-syms-internal)
ExtractSymbolsGenerateNamespaceHeader(
    ${PROJECT_NAME}-syms-internal
    ${PROJECT_NAME}-syms
    "sym"
    "private_symbols.awk"
    "${PROJECT_NAME}-syms-internal.dir"
    "private_namespace.h"
)

if(JEMALLOC_BUILD_SHARED)
    # This must be set *before* the OBJECT library is created below
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Bring sources into an object library so we can include elsewhere.
# (also so we don't compile twice for both -static and -shared)
add_library(${PROJECT_NAME} OBJECT ${JEMALLOC_CMAKE_SOURCES})
add_dependencies(${PROJECT_NAME} ${PROJECT_NAME}-syms)

if(JEMALLOC_BUILD_SHARED)
    # Now add shared library. Needed for tests and benchmarks or preloading.
    add_library(${PROJECT_NAME}-shared SHARED $<TARGET_OBJECTS:${PROJECT_NAME}>)

    if(clock_gettime_in_librt)
        target_link_libraries(${PROJECT_NAME}-shared rt)
    endif()

    if(has_libpthread)
        target_link_libraries(${PROJECT_NAME}-shared pthread)
    endif()

    if(has_libm)
        target_link_libraries(${PROJECT_NAME}-shared m)
    endif()

    if(JEMALLOC_HAVE_DLSYM)
        target_link_libraries(${PROJECT_NAME}-shared dl)
    endif()

    set_target_properties(${PROJECT_NAME}-shared PROPERTIES
        OUTPUT_NAME ${PROJECT_NAME})
endif()

if(JEMALLOC_BUILD_STATIC)
    add_library(${PROJECT_NAME}-static STATIC $<TARGET_OBJECTS:${PROJECT_NAME}>)

    # Set the name of library filename to be "libjemalloc.*"
    # instead of their full CMake libary names of libjemalloc-static.a, etc.
    set_target_properties(${PROJECT_NAME}-static PROPERTIES
        OUTPUT_NAME ${PROJECT_NAME})
endif()

if(JEMALLOC_PREFIX)
    if(JEMALLOC_BUILD_STATIC)
        target_compile_definitions(${PROJECT_NAME}-static PRIVATE JEMALLOC_MANGLE)
    endif()

    if(JEMALLOC_BUILD_SHARED)
        target_compile_definitions(${PROJECT_NAME}-shared PRIVATE JEMALLOC_MANGLE)
    endif()
endif()

# Right now we use Debug mode to imply ENABLE BUILDING ALL TESTS
# We could add an explicit flag for JEMALLOC_TESTS_ENABLE later if the Debug
# option doesn't feel like the correct place for enabling tests.
if(CMAKE_BUILD_TYPE MATCHES Debug)
    # JET setup
    add_library(${PROJECT_NAME}-jet OBJECT ${JEMALLOC_CMAKE_SOURCES})
    add_library(${PROJECT_NAME}-jet-static STATIC
        $<TARGET_OBJECTS:${PROJECT_NAME}-jet>)

    target_compile_definitions(${PROJECT_NAME}-jet
        PRIVATE
        JEMALLOC_JET
    )

    # JET sym extraction setup
    add_library(${PROJECT_NAME}-jet-syms-internal
        STATIC
        ${JEMALLOC_CMAKE_SOURCES}
    )

    add_dependencies(${PROJECT_NAME}-jet-syms-internal
        jemalloc-private-symbols-awk-jet)

    target_compile_definitions(${PROJECT_NAME}-jet-syms-internal
        PRIVATE
        JEMALLOC_JET)

    set_target_properties(${PROJECT_NAME}-jet-syms-internal PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "syms/")

    add_custom_target(${PROJECT_NAME}-jet-syms DEPENDS
        ${PROJECT_NAME}-jet-syms-internal
    )

    # Attach symbol extraction as dependency of JET
    add_dependencies(${PROJECT_NAME}-jet ${PROJECT_NAME}-jet-syms)

    ExtractSymbolsGenerateNamespaceHeader(
        ${PROJECT_NAME}-jet-syms-internal
        ${PROJECT_NAME}-jet-syms
        "sym.jet"
        "private_symbols_jet.awk"
        "${PROJECT_NAME}-jet-syms-internal.dir"
        "private_namespace_jet.h"
    )
endif()
