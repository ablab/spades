############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(hammer CXX)

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

if (SPADES_STATIC_BUILD)
  set_target_properties(hammer PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

