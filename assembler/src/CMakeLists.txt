# -*- cmake -*-
############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# Hard prerequisites
cmake_minimum_required(VERSION 2.8.12)

# allow more human readable "if then else" constructs
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE) 

# modules to find additional libraries (like zlib, boost and etc)
# Add path for custom modules
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
# Define various dirs
set(SPADES_MAIN_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SPADES_MODULES_DIR ${SPADES_MAIN_SRC_DIR}/common)
set(SPADES_MAIN_INCLUDE_DIR ${SPADES_MAIN_SRC_DIR}/include)
set(SPADES_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(SPADES_TOOLS_BINARY_DIR ${SPADES_BINARY_DIR}/bin)
set(SPADES_BUILT_INCLUDE_DIR ${SPADES_BINARY_DIR}/include)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SPADES_TOOLS_BINARY_DIR})
set(EXT_DIR "${CMAKE_SOURCE_DIR}/../ext")
set(SPADES_CFG_DIR "${CMAKE_SOURCE_DIR}/../configs")

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

message("Building SPAdes from ${SPADES_GIT_REFSPEC}, sha: ${SPADES_GIT_SHA1}")

configure_file("${SPADES_MAIN_INCLUDE_DIR}/config.hpp.in"
               "${SPADES_BUILT_INCLUDE_DIR}/config.hpp")
configure_file("${SPADES_MAIN_INCLUDE_DIR}/version.hpp.in"
               "${SPADES_BUILT_INCLUDE_DIR}/version.hpp")

# Build external dependencies (if any)
add_subdirectory("${EXT_DIR}/src" "${Project_BINARY_DIR}/ext")
add_subdirectory("${EXT_DIR}/tools" "${Project_BINARY_DIR}/ext_tools")

# Setup warnings
add_compile_options(-Wall -Wextra -Wconversion -Wno-sign-conversion -Wno-long-long -Wwrite-strings)
if (NOT OPENMP_FOUND)
  add_compile_options(-Wno-unknown-pragmas)
endif()

# sub projects
add_subdirectory(common)
add_subdirectory(projects)
add_subdirectory(spades_pipeline)


# Main pipeline script
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/../spades.py"
        DESTINATION bin
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/../dipspades.py"
        DESTINATION bin
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/../metaspades.py"
        DESTINATION bin
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/../plasmidspades.py"
        DESTINATION bin
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/../rnaspades.py"
        DESTINATION bin
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/../truspades.py"
        DESTINATION bin
        COMPONENT runtime)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../spades_init.py"
        DESTINATION bin
        COMPONENT runtime)
# FIXME: Should we just configure spades_init instead?
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../VERSION"
        DESTINATION share/spades
        COMPONENT runtime)
# Test datasets
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../test_dataset"
        DESTINATION share/spades)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../test_dataset_truspades"
        DESTINATION share/spades)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../test_dataset_plasmid"
        DESTINATION share/spades)
# manual, LICENSE, and GPLv2
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../manual.html"
        DESTINATION share/spades
        COMPONENT runtime)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../dipspades_manual.html"
        DESTINATION share/spades
        COMPONENT runtime)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../truspades_manual.html"
        DESTINATION share/spades
        COMPONENT runtime)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../rnaspades_manual.html"
        DESTINATION share/spades
        COMPONENT runtime)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE"
        DESTINATION share/spades
        COMPONENT runtime)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../GPLv2.txt"
        DESTINATION share/spades
        COMPONENT runtime)

# Random stuff
include(CMakeListsInternal.txt)

# Packaging
include(pack)

# Workaround for poor CLion users
file(GLOB_RECURSE clion_all_headers
  ${CMAKE_SOURCE_DIR}/*.hpp ${CMAKE_SOURCE_DIR}/*.h)
add_custom_target(all_clion SOURCES ${clion_all_headers})
