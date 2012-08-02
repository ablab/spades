############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

cmake_minimum_required (VERSION 2.8)
 
# allow more human readable "if then else" constructs
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE) 

# modules to find additional libraries (like zlib, boost and etc)
# Add path for custom modules
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

# Define various dirs
set(SPADES_MAIN_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SPADES_MAIN_INCLUDE_DIR ${SPADES_MAIN_SRC_DIR}/include)
set(SPADES_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(SPADES_TOOLS_BINARY_DIR ${SPADES_BINARY_DIR}/bin)
set(SPADES_BUILT_INCLUDE_DIR ${SPADES_BINARY_DIR}/include)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SPADES_TOOLS_BINARY_DIR})

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

# Check for OpenMP and add flags if necessary
find_package(OpenMP)
if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# Use parallel libstdc++ if possible
  add_definitions(-DUSE_GLIBCXX_PARALLEL=1)
else ()
  if (NOT APPLE)
    message(FATAL_ERROR "SPAdes requires OpenMP to be available")
  endif()
endif()

add_definitions(-std=c++0x)
add_definitions(-Wall -Wno-long-long -Wno-deprecated)
if (NOT OPENMP_FOUND)
  add_definitions(-Wno-unknown-pragmas)
endif()

# Use libc++ with clang due to C++11 mode
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  message("Making Debug Configuration...")

  add_definitions(-g3)
  add_definitions(-D_GLIBCXX_DEBUG)
else()
  message("Making Release Configuration...")

  if (${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
    add_definitions(-g3)
  else()
    add_definitions(-g0)
  endif()

  add_definitions(-O2)
  add_definitions(-DNDEBUG)
endif()


# common dependencies
find_package(ZLIB)

set(Boost_USE_MULTITHREADED ON)
find_package(Boost COMPONENTS system filesystem REQUIRED)

set(EXT_DIR "${CMAKE_SOURCE_DIR}/../ext")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${SPADES_MAIN_INCLUDE_DIR} ${SPADES_BUILT_INCLUDE_DIR})
include_directories("${EXT_DIR}/include")
include_directories("${ZLIB_INCLUDE_DIRS}")
include_directories("${Boost_INCLUDE_DIRS}")

set(COMMON_LIBRARIES ${ZLIB_LIBRARIES} ${Boost_LIBRARIES})

# sub projects
add_subdirectory(io)
add_subdirectory(debruijn)
add_subdirectory(hammer)

# Random stuff
include(CMakeListsInternal.txt)
