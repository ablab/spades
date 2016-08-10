############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(common_modules CXX)

add_subdirectory(pipeline)
add_subdirectory(assembly_graph)
add_subdirectory(basic/sequence)
add_subdirectory(utils/math)
add_subdirectory(modules/path_extend)
add_subdirectory(modules)
add_subdirectory(paired_info)
add_subdirectory(stages)
add_subdirectory(utils)
add_subdirectory(io)
add_subdirectory(utils/mph_index)

add_library(common_modules STATIC empty.cpp)

target_link_libraries(common_modules assembly_graph input sequence pipeline math_module path_extend paired_info stages utils mph_index modules)
