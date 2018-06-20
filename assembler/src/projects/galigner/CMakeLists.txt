############################################################################
# Copyright (c) 2015-2018 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(galigner CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(mapping_printer STATIC mapping_printer.cpp)

add_executable(longreads_aligner align_longreads.cpp)
target_link_libraries(longreads_aligner mapping_printer common_modules cityhash bwa edlib graphio ${COMMON_LIBRARIES})

add_executable(idealreads_aligner align_idealreads.cpp)
target_link_libraries(idealreads_aligner common_modules cityhash bwa edlib graphio ${COMMON_LIBRARIES})

add_executable(rna16S_aligner align_16S.cpp)
target_link_libraries(rna16S_aligner common_modules cityhash bwa edlib graphio ${COMMON_LIBRARIES})