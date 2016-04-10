############################################################################
# Copyright (c) 2015-2016 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

cmake_minimum_required(VERSION 2.8)

project(kmer_count_filter CXX)

include_directories(kmc_api)
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${EXT_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/debruijn)

add_executable(kmer_multiplicity_counter
        kmc_api/kmc_file.cpp
        kmc_api/kmer_api.cpp
        kmc_api/mmer.cpp
        kmer_multiplicity_counter.cpp)

target_link_libraries(kmer_multiplicity_counter input dev_support ${COMMON_LIBRARIES})

add_executable(annotation_propagator
               propagate.cpp)

target_link_libraries(annotation_propagator spades_modules nlopt BamTools ssw yaml-cpp ${COMMON_LIBRARIES})

add_executable(read_binning
               read_binning.cpp)

target_link_libraries(read_binning spades_modules nlopt BamTools ssw yaml-cpp ${COMMON_LIBRARIES})

add_executable(stats
               stats.cpp
               ../../modules/data_structures/sequence/genome_storage.cpp)

target_link_libraries(stats spades_modules nlopt BamTools ssw yaml-cpp ${COMMON_LIBRARIES})

add_executable(contig_abundance_counter
               contig_abundance_counter.cpp)

target_link_libraries(contig_abundance_counter spades_modules nlopt BamTools ssw yaml-cpp ${COMMON_LIBRARIES})
