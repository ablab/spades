############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(gfaio CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(gfaio STATIC
            gfa_reader.cpp)
target_link_libraries(gfaio gfa1)
