############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades_modules CXX)

add_subdirectory(pipeline)
add_subdirectory(assembly_graph/graph_support)
add_subdirectory(data_structures/sequence)
add_subdirectory(math)
add_subdirectory(algorithms/path_extend)
add_subdirectory(algorithms)
add_subdirectory(paired_info)
add_subdirectory(stages)
add_subdirectory(dev_support)
add_subdirectory(io)
add_subdirectory(data_structures/mph_index)

add_library(spades_modules STATIC)

target_link_libraries(spades_modules graph_support input sequence pipeline math_module path_extend paired_info stages dev_support mph_index algorithms)
