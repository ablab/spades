############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(truseq_analysis CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${CMAKE_SOURCE_DIR}/debruijn)
include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(truseq_analysis
               ../debruijn/kmer_coverage_model.cpp
               ../debruijn/config_struct.cpp
               ../debruijn/path_extend/pe_config_struct.cpp
               ../debruijn/stage.cpp
               ../debruijn/construction.cpp
               ../debruijn/gap_closer.cpp
               ../debruijn/simplification.cpp
               ../debruijn/mismatch_correction.cpp
               ../debruijn/pair_info_count.cpp
               ../debruijn/distance_estimation.cpp
               ../debruijn/genomic_info_filler.cpp
               ../debruijn/config_struct.cpp
          main.cpp alignment_analyser.cpp AlignmentAnalyserNew.cpp consistent_mapping.cpp analysis_pipeline.cpp)

target_link_libraries(truseq_analysis input mph_index yaml-cpp nlopt ${COMMON_LIBRARIES} ssw)

