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
    vector<EdgeId> intermediate_path;
    int return_code = -1;
};

struct MappingPoint {
    size_t seq_pos;
    size_t edge_pos;
};

struct PathRange {
    MappingPoint path_start;
    MappingPoint path_end;
};

struct GraphPosition {
    EdgeId edgeid;
    int position;

    GraphPosition(EdgeId edgeid_, int position_):
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

    void PrepareInitialState(omnigraph::MappingPath<debruijn_graph::EdgeId> &path,
                             const Sequence &s,
                             bool forward,
                             Sequence &ss,
                             EdgeId &start_e, int &start_pos, int &start_pos_seq) const;

    void UpdatePath(vector<debruijn_graph::EdgeId> &path,
                    std::vector<EdgeId> &ans,
                    int end_pos, int end_pos_seq, PathRange &range, bool forward) const;

  public:

    GapFiller(const debruijn_graph::Graph &g,
              const debruijn_graph::config::pacbio_processor &pb_config,
              const GapClosingConfig &gap_cfg):
        g_(g), pb_config_(pb_config), gap_cfg_(gap_cfg) {}

    GapFillerResult Run(const string &s,
                        const GraphPosition &start_pos,
                        const GraphPosition &end_pos,
                        int path_min_length, int path_max_length) const;

    GapFillerResult Run(omnigraph::MappingPath<debruijn_graph::EdgeId> &bwa_hits,
                        vector<debruijn_graph::EdgeId> &path,
                        const Sequence &s, bool forward, PathRange &range, int &return_code) const;

  private:
    const debruijn_graph::Graph &g_;
    const debruijn_graph::config::pacbio_processor pb_config_;
    const GapClosingConfig gap_cfg_;
};


} // namespace sensitive_aligner

