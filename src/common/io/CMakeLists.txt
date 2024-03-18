############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(input CXX)

add_subdirectory(binary)
add_subdirectory(sam)
set(IO_LIBS samio binary_io)
if (SPADES_USE_NCBISDK)
  add_subdirectory(sra)
  set(IO_LIBS ${IO_LIBS} sraio)
endif()

add_library(input STATIC
            reads/parser.cpp
            reads/paired_readers.cpp
            reads/binary_converter.cpp
            reads/binary_streams.cpp
            reads/io_helper.cpp
            dataset_support/read_converter.cpp
            dataset_support/dataset_readers.cpp)
target_link_libraries(input ${IO_LIBS} zlibstatic)

add_subdirectory(graph)
