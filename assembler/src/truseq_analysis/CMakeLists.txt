############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(truseq_analysis CXX)

add_executable(truseq_analysis
<<<<<<< HEAD
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
=======
               main.cpp
               alignment_analyser.cpp AlignmentAnalyserNew.cpp consistent_mapping.cpp analysis_pipeline.cpp)
>>>>>>> 0897e96... Introduce new 'debruijn' project which contains bulk of SPAdes core.

target_include_directories(truseq_analysis PRIVATE ${CMAKE_SOURCE_DIR}/debruijn)
target_link_libraries(truseq_analysis debruijn)

