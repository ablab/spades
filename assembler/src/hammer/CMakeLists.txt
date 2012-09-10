############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(hammer CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(hammer
   main.cpp
   hammer_tools.cpp
   hamcluster.cpp
   kmer_cluster.cpp
   kmer_data.cpp
   position_read.cpp
   globals.cpp
   config_struct_hammer.cpp
   mph_index/mph_index.cc
   mph_index/MurmurHash3.cpp
   mph_index/mph_bits.cc
   mph_index/trigraph.cc
   )

if (SPADES_BUILD_INTERNAL)
  add_subdirectory(quake_count)
  add_subdirectory(gen_test_data)
else()
  add_subdirectory(quake_count EXCLUDE_FROM_ALL)
  add_subdirectory(gen_test_data EXCLUDE_FROM_ALL)
endif()

target_link_libraries(hammer ${COMMON_LIBRARIES} input)

if (SPADES_STATIC_BUILD)
  set_target_properties(hammer PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS hammer
        RUNTIME DESTINATION bin)
