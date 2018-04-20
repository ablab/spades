############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(input CXX)

add_library(input STATIC
            reads/parser.cpp
            sam/read.cpp
            sam/sam_reader.cpp)

target_link_libraries(input BamTools samtools)

add_subdirectory(graph)
