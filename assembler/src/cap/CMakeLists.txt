############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(cap CXX)

include_directories(${CMAKE_SOURCE_DIR}/debruijn)

add_executable(cap
 ${EXT_DIR}/include/teamcity_boost/teamcity_boost.cpp
 ${EXT_DIR}/include/teamcity_boost/teamcity_messages.cpp
 tools.cpp)

#if (DEFINED static_build)
#  set_target_properties(cap PROPERTIES LINK_SEARCH_END_STATIC 1)
#endif()

target_link_libraries(cap ${COMMON_LIBRARIES} input)
