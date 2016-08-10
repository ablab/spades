############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
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
               config_struct_hammer.cpp
               read_corrector.cpp
               expander.cpp)

#  add_subdirectory(quake_count)
#  add_subdirectory(gen_test_data)

target_link_libraries(hammer input utils mph_index pipeline BamTools format ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(hammer PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS hammer
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/hammer"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info")
