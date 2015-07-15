############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(debruijn_test CXX)

include_directories(${EXT_DIR}/include/ConsensusCore)
include_directories("${CMAKE_SOURCE_DIR}/debruijn")

add_executable(debruijn_test
 ${EXT_DIR}/include/teamcity_boost/teamcity_boost.cpp
 ${EXT_DIR}/include/teamcity_boost/teamcity_messages.cpp
 ${CMAKE_SOURCE_DIR}/debruijn/config_struct.cpp
 ${CMAKE_SOURCE_DIR}/debruijn/path_extend/pe_config_struct.cpp
 ${CMAKE_SOURCE_DIR}/debruijn/path_extend/bidirectional_path.cpp
 test.cpp)

target_link_libraries(debruijn_test ${COMMON_LIBRARIES} ConsensusCore input mph_index nlopt ssw format ${ZLIB_LIBRARIES})

add_executable(component_generator
 generate_component.cpp)

target_link_libraries(component_generator ${COMMON_LIBRARIES} ConsensusCore input mph_index nlopt ssw format ${ZLIB_LIBRARIES})

add_executable(reference_fixer
 reference_fixer.cpp)
 
target_link_libraries(reference_fixer ${COMMON_LIBRARIES} ConsensusCore input mph_index nlopt ssw format ${ZLIB_LIBRARIES})
