############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(modules CXX)

add_library(modules STATIC
            chromosome_remover.cpp
            genome_consistance_checker.cpp
            path_extend/scaff_supplementary.cpp)

target_link_libraries(modules sequence bwa assembly_graph alignment)
