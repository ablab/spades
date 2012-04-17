############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(hammer CXX C)

find_package(OpenMP)

# bug in cmake, fixed in 2.6.4 (patch 8)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp")
# bug, problem occured in Ubuntu 12.04
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_definitions(${OpenMP_CXX_FLAGS})
add_definitions(-fPIC)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(hammer
   main.cpp
   hammer_tools.cpp
   kmer_cluster.cpp
   position_read.cpp
   globals.cpp
   subkmers.cpp
   kmerno.cpp
   kmer_stat.cpp
   config_struct_hammer.cpp
   read_reconstruction.cpp )
   
add_subdirectory(quake_count)
#add_subdirectory(quake_enhanced)
add_subdirectory(quake_correct)
add_subdirectory(gen_test_data)

target_link_libraries(hammer ${COMMON_LIBRARIES} input)

if (DEFINED static_build)
  set_target_properties(hammer PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

