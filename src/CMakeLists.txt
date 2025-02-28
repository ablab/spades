# -*- cmake -*-
############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# Hard prerequisites
cmake_minimum_required(VERSION 3.16)

project(SpadesToolkit
        VERSION 3.16.0
        LANGUAGES C CXX)

# allow more human readable "if then else" constructs
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# modules to find additional libraries (like zlib, boost and etc)
# Add path for custom modules
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
# Add some utility macro
include(utils)

# Define various dirs
set(SPADES_MAIN_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SPADES_MAIN_PROJ_DIR ${SPADES_MAIN_SRC_DIR}/projects)
set(SPADES_MODULES_DIR ${SPADES_MAIN_SRC_DIR}/common)
set(SPADES_MAIN_INCLUDE_DIR ${SPADES_MAIN_SRC_DIR}/include)
set(SPADES_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(SPADES_TOOLS_BINARY_DIR ${SPADES_BINARY_DIR}/bin)
set(SPADES_BUILT_INCLUDE_DIR ${SPADES_BINARY_DIR}/include)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SPADES_TOOLS_BINARY_DIR})
set(EXT_DIR "${CMAKE_SOURCE_DIR}/../ext")

# Uncomment for gprof profiling
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
#SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

# Everything option-dependent
include(options)

# Sanity check our source directory to make sure that we are not trying to
# generate an in-tree build
if ((CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR) OR
    (CMAKE_SOURCE_DIR STREQUAL "${CMAKE_BINARY_DIR}/src"))
  message(FATAL_ERROR "In-source builds are not allowed.
Please create a directory and run cmake from there, passing the path
to this source directory as the last argument.
This process created the file `CMakeCache.txt' and the directory `CMakeFiles'.
Please delete them.")
endif()

# Handle commond deps
include(deps)

# Add flags if necessary
include(flags)

# Handle common includes
include(includes)

# Handle common libs
include(libs)

include(GetGitRevisionDescription)
get_git_head_revision(SPADES_GIT_REFSPEC SPADES_GIT_SHA1)
if ((NOT SPADES_GIT_REFSPEC) OR (SPADES_GIT_REFSPEC STREQUAL "GITDIR-NOTFOUND"))
  set(SPADES_GIT_REFSPEC "N/A")
endif()
if (NOT SPADES_GIT_SHA1)
  set(SPADES_GIT_SHA1 "N/A")
endif()

file(STRINGS "${CMAKE_SOURCE_DIR}/../VERSION" VERSION_STRING LIMIT_COUNT 1)
string(REGEX MATCH "([a-zA-Z]*)-\?([0-9]+)\.([0-9]+)\.([0-9]+)(.*)" ver "${VERSION_STRING}")
if (ver STREQUAL "")
  message(FATAL_ERROR "Error to parse version string from ${VERSION_STRING}")
endif()
if (CMAKE_MATCH_1 STREQUAL "")
  set(SPADES_FLAVOUR "SPAdes")
else()
  set(SPADES_FLAVOUR "${CMAKE_MATCH_1}")
endif()
set(SPADES_VERSION_MAJOR "${CMAKE_MATCH_2}")
set(SPADES_VERSION_MINOR "${CMAKE_MATCH_3}")
set(SPADES_VERSION_PATCH "${CMAKE_MATCH_4}")
set(SPADES_VERSION_SUFFIX "${CMAKE_MATCH_5}")
set(SPADES_PACKAGE_VERSION "${SPADES_VERSION_MAJOR}.${SPADES_VERSION_MINOR}.${SPADES_VERSION_PATCH}${SPADES_VERSION_SUFFIX}")
message(STATUS "Building ${SPADES_FLAVOUR} version ${SPADES_PACKAGE_VERSION}")

# Testing
include(test)

message("Building SPAdes from ${SPADES_GIT_REFSPEC}, sha: ${SPADES_GIT_SHA1}")

configure_file("${SPADES_MAIN_INCLUDE_DIR}/config.hpp.in"
               "${SPADES_BUILT_INCLUDE_DIR}/config.hpp")
configure_file("${SPADES_MAIN_INCLUDE_DIR}/version-gen.hpp.in"
               "${SPADES_BUILT_INCLUDE_DIR}/version-gen.hpp")

# Build external dependencies (if any)
add_subdirectory("${EXT_DIR}/src" "${PROJECT_BINARY_DIR}/ext")
add_subdirectory("${EXT_DIR}/include" "${PROJECT_BINARY_DIR}/ext-include-only")

# Setup warnings
add_compile_options(-Wall -Wextra -Wconversion -Wno-sign-conversion -Wno-long-long -Wwrite-strings)
if (NOT OPENMP_FOUND)
  add_compile_options(-Wno-unknown-pragmas)
endif()

# Looking for README.md
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../README.md")
  #installing from archive
  set(README_MD_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../README.md")
else()
  #installing from the repo
  set(README_MD_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../README.md")
endif()

# sub projects
include(proj)
add_subdirectory(common)
add_subdirectory(projects)

# FIXME: Should we just configure spades_init instead?
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../VERSION"
        DESTINATION share/spades
        COMPONENT common)
# manual, LICENSE, and GPLv2
install(FILES "${README_MD_PATH}"
        DESTINATION share/spades
        COMPONENT spades-docs)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE"
        DESTINATION share/spades
        COMPONENT common)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../GPLv2.txt"
        DESTINATION share/spades
        COMPONENT common)

# Random stuff
include(CMakeListsInternal.txt)

# Packaging
include(pack)

# Workaround for poor CLion users
file(GLOB_RECURSE clion_all_headers
  ${CMAKE_SOURCE_DIR}/*.hpp ${CMAKE_SOURCE_DIR}/*.h)
add_custom_target(all_clion SOURCES ${clion_all_headers})
