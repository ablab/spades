############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(truseq_analysis CXX)

add_executable(truseq_analysis
               main.cpp
               alignment_analyser.cpp AlignmentAnalyserNew.cpp consistent_mapping.cpp analysis_pipeline.cpp)

target_link_libraries(truseq_analysis common_modules ${COMMON_LIBRARIES})

