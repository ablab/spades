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
target_link_libraries(debruijn_test common_modules cityhash ssw edlib ${COMMON_LIBRARIES})

add_executable(component_generator generate_component.cpp)
target_link_libraries(component_generator common_modules cityhash ${COMMON_LIBRARIES})

add_executable(reference_fixer reference_fixer.cpp)
target_link_libraries(reference_fixer common_modules cityhash ${COMMON_LIBRARIES})

add_executable(short_edge_contractor short_edge_contractor.cpp)
target_link_libraries(short_edge_contractor common_modules cityhash ${COMMON_LIBRARIES})

add_executable(sequence_threader thread_sequences.cpp)
target_link_libraries(sequence_threader graphio common_modules ${COMMON_LIBRARIES})
