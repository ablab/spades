############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(pipeline CXX)

add_library(pipeline STATIC
            graph_pack.cpp
            graph_pack_helpers.cpp
            sequence_mapper_gp_api.cpp
            stage.cpp)

target_link_libraries(pipeline binary_io path_extend input llvm-support library configs alignment)