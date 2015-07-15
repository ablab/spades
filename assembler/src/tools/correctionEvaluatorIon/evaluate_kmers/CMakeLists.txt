############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(kmer_evaluator CXX)

include_directories(../../../hammer-it)

add_executable(kmer_evaluator
               ../../../hammer-it/kmer_data.cpp
               ../../../hammer-it/err_helper_table.cpp
               kmer_evaluator.cpp)

target_link_libraries(kmer_evaluator input mph_index yaml-cpp ${COMMON_LIBRARIES})
