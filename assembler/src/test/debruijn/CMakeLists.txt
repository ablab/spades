############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(debruijn_test CXX)

include_directories("${CMAKE_SOURCE_DIR}/debruijn")

add_executable(debruijn_test
 ${EXT_DIR}/include/teamcity_boost/teamcity_boost.cpp
 ${EXT_DIR}/include/teamcity_boost/teamcity_messages.cpp
 ${CMAKE_SOURCE_DIR}/debruijn/kmer_coverage_model.cpp
 test.cpp)

target_link_libraries(debruijn_test ${COMMON_LIBRARIES} input mph_index otkpp)
