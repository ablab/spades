############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(debruijn_test CXX)

add_executable(debruijn_test
               ${EXT_DIR}/include/teamcity_boost/teamcity_boost.cpp
               ${EXT_DIR}/include/teamcity_boost/teamcity_messages.cpp
               test.cpp)
target_include_directories(debruijn_test PRIVATE "${CMAKE_SOURCE_DIR}/debruijn")
target_link_libraries(debruijn_test debruijn)

add_executable(component_generator generate_component.cpp)
target_link_libraries(component_generator debruijn)
target_include_directories(component_generator PRIVATE "${CMAKE_SOURCE_DIR}/debruijn")

add_executable(reference_fixer reference_fixer.cpp)
target_link_libraries(reference_fixer debruijn)
target_include_directories(reference_fixer PRIVATE "${CMAKE_SOURCE_DIR}/debruijn")
