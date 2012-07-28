############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(hammer CXX C)

find_package(OpenMP)

if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# Use parallel libstdc++ if possible
  add_definitions(-DUSE_GLIBCXX_PARALLEL=1)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(hammer
   main.cpp
   half.cpp
   hammer_tools.cpp
   hammer_stats.cpp
   hamcluster.cpp
   kmer_cluster.cpp
   position_read.cpp
   globals.cpp
   kmerno.cpp
   kmer_stat.cpp
   config_struct_hammer.cpp
   )
   
add_subdirectory(quake_count)
#add_subdirectory(quake_enhanced)
#add_subdirectory(quake_correct)
add_subdirectory(gen_test_data)

target_link_libraries(hammer ${COMMON_LIBRARIES} input)

if (DEFINED static_build)
  set_target_properties(hammer PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

