############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(debruijn_test CXX)

include_directories("${CMAKE_SOURCE_DIR}/debruijn")

add_executable(debruijn_test
 ${EXT_DIR}/include/teamcity_boost/teamcity_boost.cpp
 ${EXT_DIR}/include/teamcity_boost/teamcity_messages.cpp
 ${CMAKE_SOURCE_DIR}/debruijn/config_struct.cpp
 ${CMAKE_SOURCE_DIR}/debruijn/path_extend/pe_config_struct.cpp
 test.cpp)
 
include_directories(${EXT_DIR}/include/ConsensusCore)
target_link_libraries(debruijn_test ${COMMON_LIBRARIES} ConsensusCore input mph_index otkpp ssw ${ZLIB_LIBRARIES})
