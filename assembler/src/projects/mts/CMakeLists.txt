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

target_link_libraries(kmer_multiplicity_counter common_modules utils input getopt_pp ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(kmer_multiplicity_counter PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

add_executable(prop_binning
               annotation.cpp
               propagate.cpp
               read_binning.cpp
               prop_binning.cpp
               gzstream/gzstream.C)

target_link_libraries(prop_binning common_modules nlopt BamTools ssw getopt_pp ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(prop_binning PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

add_executable(stats
               annotation.cpp
               stats.cpp)

target_link_libraries(stats common_modules nlopt BamTools ssw getopt_pp ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(stats PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

add_executable(contig_abundance_counter
               contig_abundance_counter.cpp
               contig_abundance.cpp)

target_link_libraries(contig_abundance_counter common_modules nlopt BamTools ssw getopt_pp ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(contig_abundance_counter PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()
