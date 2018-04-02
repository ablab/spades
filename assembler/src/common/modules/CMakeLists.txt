############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(modules CXX)

add_library(modules STATIC
            genome_consistance_checker.cpp
            alignment/gap_info.cpp
            alignment/bwa_index.cpp)

target_link_libraries(modules bwa)
