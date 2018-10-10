#pragma once

#include "pipeline/configs/aligner_config.hpp"
#include "assembly_graph/paths/path_utils.hpp"
#include "assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "assembly_graph/paths/mapping_path.hpp"
#include "assembly_graph/paths/path_processor.hpp"

#include "modules/alignment/pacbio/gap_dijkstra.hpp"

namespace sensitive_aligner {
using debruijn_graph::EdgeId;
using debruijn_graph::VertexId;

struct GapFillerResult {
    int score = std::numeric_limits<int>::max();
    vector<EdgeId> full_intermediate_path;
    int return_code = -1;
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

    GraphPosition(EdgeId edgeid_, size_t position_):
        edgeid(edgeid_), position(position_) {}
};

class GapFiller {
    std::string PathToString(const vector<EdgeId>& path) const;

    GapFillerResult BestScoredPathDijkstra(const string &s,
                                           const GraphPosition &start_pos,
                                           const GraphPosition &end_pos,
                                           int path_max_length, int score) const;

    GapFillerResult BestScoredPathBruteForce(const string &seq_string,
            const GraphPosition &start_pos,
            const GraphPosition &end_pos,
            int path_min_length, int path_max_length) const;

    void Revert(Sequence &ss, GraphPosition &start_pos) const;

    void UpdatePath(vector<debruijn_graph::EdgeId> &path,
                    std::vector<EdgeId> &ans,
                    MappingPoint p, PathRange &range, bool forward, GraphPosition &old_start_pos) const;

  public:

    GapFiller(const debruijn_graph::Graph &g,
              const debruijn_graph::config::pacbio_processor &pb_config,
              const GapClosingConfig &gap_cfg,
              const EndsClosingConfig &ends_cfg):
        g_(g), pb_config_(pb_config), gap_cfg_(gap_cfg), ends_cfg_(ends_cfg) {}

    GapFillerResult Run(const string &s,
                        const GraphPosition &start_pos,
                        const GraphPosition &end_pos,
                        int path_min_length, int path_max_length) const;

    GapFillerResult Run(Sequence &s,
                        GraphPosition &start_pos,
                        bool forward, 
                        vector<debruijn_graph::EdgeId> &path, 
                        PathRange &range) const;

  private:
    const debruijn_graph::Graph &g_;
    const debruijn_graph::config::pacbio_processor pb_config_;
    const GapClosingConfig gap_cfg_;
    const EndsClosingConfig ends_cfg_;
};


} // namespace sensitive_aligner

