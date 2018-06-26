#pragma once

namespace rna16S_mapping {

struct RnaAlignerConfig {
    // general
    int K;
    string path_to_graphfile;
    string path_to_sequences;
    string path_to_primers;

    double path_limit_stretching;
    int gap_min_ed;
    int res_min_ed;
    int position_diff;
    int alignment_min_length;
    int primer_max_ed;
    int graph_max_ed;

    graph_aligner::GapClosingConfig gap_cfg;
};

}