//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#pragma once

#include <cstdlib>

namespace debruijn_graph {
namespace config {

struct pacbio_processor {
    size_t bwa_length_cutoff      = 200;
    double compression_cutoff     = 0.6;
    double path_limit_stretching  = 1.3;
    double path_limit_pressing    = 0.7;
    size_t max_path_in_dijkstra   = 15000;
    size_t max_vertex_in_dijkstra = 2000;
    // gap closer
    size_t long_seq_limit           = 400;
    size_t pacbio_min_gap_quantity  = 2;
    size_t contigs_min_gap_quantity = 1;
    size_t max_contigs_gap_length   = 10000;
};

struct bwa_aligner {
    size_t min_contig_len;
};

}
}
