############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(cap CXX)

# Define minimum and maximum K
set(SPADES_MIN_K 1 CACHE INTEGER "Minimum k-mer length")
set(SPADES_MAX_K 300 CACHE INTEGER "Maximum k-mer length")
configure_file("${SPADES_MAIN_INCLUDE_DIR}/k_range.hpp.in"
               "${SPADES_BUILT_INCLUDE_DIR}/k_range.hpp")
include_directories(${CMAKE_SOURCE_DIR}/debruijn)

add_executable(cap
 ${EXT_DIR}/include/teamcity_boost/teamcity_boost.cpp
 ${EXT_DIR}/include/teamcity_boost/teamcity_messages.cpp
 tools.cpp)

#if (DEFINED static_build)
#  set_target_properties(cap PROPERTIES LINK_SEARCH_END_STATIC 1)
#endif()

target_link_libraries(cap ${COMMON_LIBRARIES} input mph_index)
