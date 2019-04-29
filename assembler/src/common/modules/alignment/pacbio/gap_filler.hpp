//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/configs/aligner_config.hpp"
#include "assembly_graph/paths/path_utils.hpp"
#include "assembly_graph/paths/mapping_path.hpp"
#include "assembly_graph/paths/path_processor.hpp"

#include "modules/alignment/bwa_index.hpp"
#include "modules/alignment/pacbio/gap_dijkstra.hpp"

namespace sensitive_aligner {
using debruijn_graph::EdgeId;
using debruijn_graph::VertexId;

struct GAlignerConfig {
    // general
    int K = -1;
    std::string path_to_graphfile = "";
    std::string path_to_sequences = "";
    alignment::BWAIndex::AlignmentMode data_type = alignment::BWAIndex::AlignmentMode::Default; // pacbio, nanopore, 16S
    std::string output_format = "tsv"; // default: tsv
    bool restore_ends = false;

    //path construction
    debruijn_graph::config::pacbio_processor pb;
    GapClosingConfig gap_cfg;
    EndsClosingConfig ends_cfg;

    GAlignerConfig(const debruijn_graph::config::pacbio_processor &pb_,
                   const alignment::BWAIndex::AlignmentMode &data_type_,
                   const GapClosingConfig &gap_cfg_ = GapClosingConfig(),
                   const EndsClosingConfig &ends_cfg_ = EndsClosingConfig())
    :data_type(data_type_),
     pb(pb_),
     gap_cfg(gap_cfg_),
     ends_cfg(ends_cfg_)
    {}

    GAlignerConfig(){}
};


struct GapFillerResult {
    int score = std::numeric_limits<int>::max();
    std::vector<EdgeId> full_intermediate_path;
    DijkstraReturnCode return_code;
};

struct MappingPoint {
    size_t seq_pos;
    size_t edge_pos;

    MappingPoint(): seq_pos(0), edge_pos(0) {}

    MappingPoint(size_t seq_pos_, size_t edge_pos_):
        seq_pos(seq_pos_), edge_pos(edge_pos_) {}

};

struct PathRange {
    MappingPoint path_start;
    MappingPoint path_end;

    PathRange() {}

    PathRange(MappingPoint a, MappingPoint b):
        path_start(a), path_end(b) {}
};

struct GraphPosition {
    EdgeId edgeid;
    size_t position;

    GraphPosition(EdgeId edgeid_ = EdgeId(), size_t position_ = -1):
        edgeid(edgeid_), position(position_) {}
};

class GapFiller {
  public:

    GapFiller(const debruijn_graph::Graph &g,
              const GAlignerConfig &cfg):
        g_(g), cfg_(cfg) {}

    GapFillerResult Run(const std::string &s,
                        const GraphPosition &start_pos,
                        const GraphPosition &end_pos,
                        int path_min_length, int path_max_length) const;

    GapFillerResult Run(Sequence &s,
                        GraphPosition &start_pos,
                        bool forward,
                        std::vector<debruijn_graph::EdgeId> &path,
                        PathRange &range) const;

  private:

    std::string PathToString(std::vector<EdgeId> &path) const;

    GapFillerResult BestScoredPathDijkstra(const std::string &s,
                                           const GraphPosition &start_pos,
                                           const GraphPosition &end_pos,
                                           int path_max_length, int score) const;

    GapFillerResult BestScoredPathBruteForce(const std::string &seq_string,
            const GraphPosition &start_pos,
            const GraphPosition &end_pos,
            int path_min_length, int path_max_length) const;

    GraphPosition ConjugatePosition(const GraphPosition &start_pos) const;

    void UpdatePath(std::vector<debruijn_graph::EdgeId> &path,
                    std::vector<EdgeId> &ans,
                    MappingPoint p, PathRange &range, bool forward, GraphPosition &old_start_pos) const;

    const debruijn_graph::Graph &g_;
    const GAlignerConfig &cfg_;
};


} // namespace sensitive_aligner

