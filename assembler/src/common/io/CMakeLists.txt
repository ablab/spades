############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(input CXX)

add_subdirectory(binary)

add_library(input STATIC
            reads/parser.cpp
            reads/paired_readers.cpp
            reads/binary_converter.cpp
            reads/binary_streams.cpp
            reads/io_helper.cpp
            dataset_support/read_converter.cpp
            dataset_support/dataset_readers.cpp
            sam/read.cpp
            sam/sam_reader.cpp)

include_directories(SYSTEM "${ZLIB_INCLUDE_DIRS}")
target_link_libraries(input binary_io BamTools samtools ${ZLIB_LIBRARIES})

add_subdirectory(graph)
