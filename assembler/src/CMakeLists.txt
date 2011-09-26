cmake_minimum_required (VERSION 2.6)
 
# allow more human readable "if then else" constructs
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE) 

# modules to find additional libraries (like log4cxx, apr and etc)
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

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

add_definitions(-g3)
add_definitions(-std=c++0x)
add_definitions(-Wall -Wno-long-long -Wno-deprecated)
add_definitions(-rdynamic)
#add_definitions(-Wall -Wconversion -Wno-long-long -Wno-deprecated)
add_definitions(-Wall -Wno-long-long -Wno-deprecated)

if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  message("Making Debug Configuration...")
  add_definitions(-D_GLIBCXX_DEBUG)
else()
  message("Making Release Configuration...")
  add_definitions(-o2)
  
  #temporarily turned off to run our assembler on large datasets
  #add_definitions(-DNDEBUG)
endif()


# common dependencies
find_package(ZLIB)
find_package(log4cxx)
find_package(statgen)
find_package(staden)
find_package(ssl)
find_package(bam)
find_package(Boost COMPONENTS signals system filesystem REQUIRED)
find_package(fftw3)

set(EXT_DIR "${CMAKE_SOURCE_DIR}/../ext")

include_directories("${CMAKE_SOURCE_DIR}/include")

include_directories(${ZLIB_INCLUDE_DIRS})
include_directories(${log4cxx_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${statgen_INCLUDE_DIRS})
include_directories(${staden_INCLUDE_DIRS})
include_directories(${ssl_INCLUDE_DIRS})
include_directories(${bam_INCLUDE_DIRS})
include_directories(${fftw3_INCLUDE_DIRS})

set(COMMON_LIBRARIES 
   ${ZLIB_LIBRARIES} 
   ${log4cxx_LIBRARIES} 
   ${Boost_LIBRARIES} 
   ${statgen_LIBRARIES} 
   ${staden_LIBRARIES} 
   ${ssl_LIBRARIES} 
   ${bam_LIBRARIES}
   ${fftw3_LIBRARIES}
   )

include_directories(${EXT_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include)

# sub projects
add_subdirectory(io)
add_subdirectory(debruijn)
add_subdirectory(abruijn)
add_subdirectory(bayes)
add_subdirectory(hammer)
add_subdirectory(common_test)
add_subdirectory(debruijn/long_contigs)
add_subdirectory(debruijn/long_contigs/utils)
