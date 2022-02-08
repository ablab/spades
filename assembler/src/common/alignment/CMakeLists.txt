############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(alignment CXX)

add_library(alignment STATIC
            gap_info.cpp
            bwa_index.cpp
            long_read_mapper.cpp
            sequence_mapper_notifier.cpp
            pacbio/gap_filler.cpp
            pacbio/gap_dijkstra.cpp
            pacbio/g_aligner.cpp
            pacbio/g_aligner.cpp)

target_link_libraries(alignment sequence bwa assembly_graph)
