############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(debruijn_examples CXX)

add_executable(component_generator generate_component.cpp)
target_link_libraries(component_generator common_modules ${COMMON_LIBRARIES})

add_executable(reference_fixer reference_fixer.cpp)
target_link_libraries(reference_fixer common_modules ${COMMON_LIBRARIES})

add_executable(short_edge_contractor short_edge_contractor.cpp)
target_link_libraries(short_edge_contractor graphio common_modules ${COMMON_LIBRARIES})

add_executable(sequence_threader thread_sequences.cpp)
target_link_libraries(sequence_threader toolchain graphio common_modules ${COMMON_LIBRARIES})
