############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(graphio CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(graphio STATIC
            gfa_reader.cpp gfa_writer.cpp
            fastg_writer.cpp)
target_link_libraries(graphio gfa1)
