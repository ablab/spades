############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(input CXX)

add_library(input STATIC
            parser.cpp
            path_helper.cpp
            copy_file.cpp
            library.cpp
            logger_impl.cpp
            sam/read.cpp
            sam/sam_reader.cpp)

target_link_libraries(input BamTools samtools yaml-cpp)

