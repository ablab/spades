############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set (VERBOSE 1)
cmake_minimum_required (VERSION 2.6)
 
# allow more human readable "if then else" constructs
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE) 

# modules to find additional libraries (like zlib, boost and etc)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

# configuration routine
# 'cfg' is short alias for CMAKE_BUILD_TYPE

if (DEFINED cfg) 
  set(CMAKE_BUILD_TYPE "${cfg}")
else() 
  set(CMAKE_BUILD_TYPE "RelWithDebInfo") 
endif()

if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug"          AND 
    NOT ${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
    
  set(CMAKE_BUILD_TYPE "RelWithDebInfo") 
endif()

if (DEFINED static_build)
  # it'll make cmake to find libraries archives, not dynamic link
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a) 
  
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")

  set(Boost_USE_STATIC_LIBS        ON)
  set(Boost_USE_STATIC_RUNTIME     ON)
endif()

find_package(OpenMP)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp")

add_definitions(${OpenMP_CXX_FLAGS})
add_definitions(-g3)
<<<<<<< HEAD
add_definitions(-static)
add_definitions(-fopenmp)
=======
>>>>>>> works a little
add_definitions(-std=c++0x)
add_definitions(-Wall -Wno-long-long -Wno-deprecated)
#add_definitions(-Wconversion)
#add_definitions(-Wall -Wno-long-long -Wno-deprecated)

if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  message("Making Debug Configuration...")
  add_definitions(-D_GLIBCXX_DEBUG)
else()
  message("Making Release Configuration...")
  add_definitions(-o2)
  
  add_definitions(-DNDEBUG)
endif()


# common dependencies
find_package(ZLIB)

set(Boost_USE_MULTITHREADED ON)
find_package(Boost COMPONENTS system filesystem iostreams serialization REQUIRED)

set(EXT_DIR "${CMAKE_SOURCE_DIR}/../ext")

include_directories("${CMAKE_SOURCE_DIR}/include")

include_directories(${ZLIB_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

set(COMMON_LIBRARIES 
   ${ZLIB_LIBRARIES} 
   ${Boost_LIBRARIES} 
   )

include_directories(${EXT_DIR}/include)
include_directories(${EXT_DIR}/include/google/)
include_directories(${CMAKE_SOURCE_DIR}/include)

# sub projects
add_subdirectory(io)
add_subdirectory(debruijn)
add_subdirectory(hammer)

include(CMakeListsInternal.txt)
<<<<<<< HEAD
<<<<<<< HEAD

=======
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread -lrt ") #-lgomp
=======
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread -lrt")
>>>>>>> works a little
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
>>>>>>> compiles
