############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(coverage_model CXX)

add_library(coverage_model STATIC
        kmer_coverage_model.cpp)

target_link_libraries(coverage_model nlopt)

