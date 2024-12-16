############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2018-2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(binary_io CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(binary_io STATIC
            graph_pack.cpp genomic_info.cpp
            )
target_link_libraries(binary_io llvm-support)
