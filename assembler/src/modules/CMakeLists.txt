############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades_modules CXX)

add_subdirectory(io_impl)
add_subdirectory(pipeline)
add_subdirectory(graph_support)
add_subdirectory(sequence)
add_subdirectory(math)
add_subdirectory(path_extend)
add_subdirectory(paired_info)
add_subdirectory(shared_stages)

add_library(spades_modules STATIC)

target_link_libraries(spades_modules graph_support input sequence pipeline math_module path_extend paired_info shared_stages cityhash)
