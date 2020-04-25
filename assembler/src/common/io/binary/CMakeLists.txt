############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(binary CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(binary STATIC
            graph_pack.cpp genomic_info.cpp
            )
