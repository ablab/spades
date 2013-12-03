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
 ${CMAKE_SOURCE_DIR}/debruijn/kmer_coverage_model.cpp
 ${CMAKE_SOURCE_DIR}/debruijn/gap_closer.cpp
 ${CMAKE_SOURCE_DIR}/debruijn/simplification.cpp
 ${CMAKE_SOURCE_DIR}/debruijn/construction.cpp
 ${CMAKE_SOURCE_DIR}/debruijn/stage.cpp
 ${CMAKE_SOURCE_DIR}/debruijn/mismatch_correction.cpp
 ${CMAKE_SOURCE_DIR}/debruijn/pacbio_aligning.cpp 
 ${CMAKE_SOURCE_DIR}/debruijn/genomic_info_filler.cpp
 ${CMAKE_SOURCE_DIR}/debruijn/pair_info_count.cpp
 ${CMAKE_SOURCE_DIR}/debruijn/distance_estimation.cpp
 ${CMAKE_SOURCE_DIR}/debruijn/repeat_resolving.cpp
 ${CMAKE_SOURCE_DIR}/debruijn/config_struct.cpp
 ${CMAKE_SOURCE_DIR}/debruijn/path_extend/pe_config_struct.cpp
 test.cpp)
 
include_directories(${EXT_DIR}/include/ConsensusCore)
target_link_libraries(debruijn_test ${COMMON_LIBRARIES} ConsensusCore input mph_index otkpp ssw ${ZLIB_LIBRARIES})
