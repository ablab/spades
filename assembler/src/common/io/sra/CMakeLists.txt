############################################################################
# Copyright (c) 2021 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(samhio CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(sraio STATIC
            sra_parser.cpp)
target_link_libraries(sraio ncbi-sdk)
