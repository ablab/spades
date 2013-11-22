############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(hammer-it CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(hammer-it
               kmer_data.cpp
               hamcluster.cpp
               subcluster.cpp
               err_helper_table.cpp
               config_struct.cpp
               expander.cpp
               main.cpp)

target_link_libraries(hammer-it input mph_index yaml-cpp ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(hammer-it PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS hammer-it
        RUNTIME DESTINATION bin)
