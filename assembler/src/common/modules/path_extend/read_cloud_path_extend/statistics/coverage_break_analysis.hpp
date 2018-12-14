#pragma once

#include "common/assembly_graph/graph_support/scaff_supplementary.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/reference_path_index.hpp"

namespace path_extend {

class GraphBreakAnalyzer {
    struct GraphBreak {
      EdgeId first_edge_;
      EdgeId second_edge_;
      VertexId first_end_;
      VertexId second_start_;

      GraphBreak(const EdgeId &first_edge_,
                 const EdgeId &second_edge_,
                 const VertexId &first_end_,
                 const VertexId &second_start_);
    };

    typedef validation::EdgeWithMapping EdgeWithMapping;
    typedef std::vector<EdgeWithMapping> MappedPath;
    typedef std::vector<EdgeId> SimplePath;
    typedef std::unordered_map<transitions::Transition, std::vector<GraphBreak>> TransitionToBreaks;

    const debruijn_graph::conj_graph_pack &gp_;

 public:
    GraphBreakAnalyzer(const conj_graph_pack &gp_);

    TransitionToBreaks GetGraphBreaks(const string &path_to_reference, size_t length_threshold) const;

    void PrintGraphBreaks(const TransitionToBreaks &transition_to_breaks, const string &output_path) const;

 private:
    std::vector<MappedPath> GetRawPaths(const string &path_to_reference) const;

    std::vector<GraphBreak> GetGraphBreaks(const SimplePath &path) const;

    std::vector<MappedPath> FixLongEdgeAlignments(const std::vector<MappedPath> &raw_paths,
                                                  const std::unordered_set<EdgeId> &long_edges,
                                                  double min_mapped_ratio) const;

    SimplePath GetBarcodedSubpath(const EdgeId &first, const EdgeId &second,
                                  const std::vector<GraphBreakAnalyzer::MappedPath> &raw_paths,
                                  const validation::ReferencePathIndex &ref_path_index, size_t min_intersection) const;
};

}