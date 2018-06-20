#pragma once

#include "pipeline/configs/aligner_config.hpp"
#include "assembly_graph/paths/path_utils.hpp"
#include "assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "assembly_graph/paths/mapping_path.hpp"
#include "assembly_graph/paths/path_processor.hpp"

#include "modules/alignment/pacbio/gap_dijkstra.hpp"

namespace gap_filler {
using debruijn_graph::EdgeId;
using debruijn_graph::VertexId;

struct GapFillerResult {
    int score = std::numeric_limits<int>::max();
    vector<EdgeId> intermediate_path;
    int return_code = -1;
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
                                           const gap_filler::GraphPosition &start_pos,
                                           const gap_filler::GraphPosition &end_pos,
                                           int path_max_length, int score) const;

    GapFillerResult BestScoredPathBruteForce(const string &seq_string,
            const gap_filler::GraphPosition &start_pos,
            const gap_filler::GraphPosition &end_pos,
            int path_min_length, int path_max_length) const;

public:

    GapFiller(const debruijn_graph::Graph &g,
              const debruijn_graph::config::pacbio_processor &pb_config,
              const gap_dijkstra::GapClosingConfig &gap_cfg):
        g_(g), pb_config_(pb_config), gap_cfg_(gap_cfg) {}

    GapFillerResult Run(const string &s,
                        const gap_filler::GraphPosition &start_pos,
                        const gap_filler::GraphPosition &end_pos,
                        int path_min_length, int path_max_length);

private:
    const debruijn_graph::Graph &g_;
    const debruijn_graph::config::pacbio_processor &pb_config_;
    const gap_dijkstra::GapClosingConfig &gap_cfg_;
};

} // namespace gap_filler

