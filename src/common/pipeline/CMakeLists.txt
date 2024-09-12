############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(pipeline CXX)

set(pipeline_src
  graph_pack.cpp
  graph_pack_helpers.cpp
  sequence_mapper_gp_api.cpp
  stage.cpp)

add_library(pipeline STATIC ${pipeline_src})
target_link_libraries(pipeline binary_io path_extend input llvm-support library configs alignment)
