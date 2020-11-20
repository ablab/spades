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
            alignment/gap_info.cpp
            alignment/bwa_index.cpp
            alignment/long_read_mapper.cpp
            alignment/sequence_mapper.cpp
            alignment/sequence_mapper_notifier.cpp
            alignment/pacbio/gap_filler.cpp
            alignment/pacbio/gap_dijkstra.cpp 
            alignment/pacbio/g_aligner.cpp 
            alignment/pacbio/g_aligner.cpp)

target_link_libraries(modules sequence bwa assembly_graph)
