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
    size_t internal_length_cutoff = 200;
    double compression_cutoff     = 0.6;
    double path_limit_stretching  = 1.3;
    double path_limit_pressing    = 0.7;
    size_t max_path_in_dijkstra   = 15000;
    size_t max_vertex_in_dijkstra = 2000;
    bool rna_filtering            = false;

    // gap closer
    size_t long_seq_limit           = 400;
    bool enable_gap_closing         = true;
    bool enable_fl_gap_closing      = true;
    size_t pacbio_min_gap_quantity  = 2;
    size_t contigs_min_gap_quantity = 1;
    size_t max_contigs_gap_length   = 10000;

    // spoa settings
    int8_t match                    = 5;
    int8_t mismatch                 = -4;
    int8_t gap_open                 = -8;
    int8_t gap_extend               = -6;
    int8_t gap_open_second          = -10;
    int8_t gap_extend_second        = -4;
};

struct bwa_aligner {
    size_t min_contig_len;
};

}
}
