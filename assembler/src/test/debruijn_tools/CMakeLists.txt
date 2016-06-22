############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(debruijn_tools CXX)

include_directories("${CMAKE_SOURCE_DIR}/debruijn")

add_executable(debruijn_tools
 ${EXT_DIR}/include/teamcity_boost/teamcity_boost.cpp
 ${EXT_DIR}/include/teamcity_boost/teamcity_messages.cpp
 tools.cpp)

target_link_libraries(debruijn_tools ConsensusCore input cityhash nlopt BamTools format ${COMMON_LIBRARIES})
