#pragma once

#include <ostream>
#include "common/modules/path_extend/read_cloud_path_extend/validation/transition_subgraph_validation.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/path_cluster_helper.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/simple_graph.hpp"
#include "common/barcode_index/cluster_storage_extractor.hpp"

namespace path_extend {

struct SubgraphInfo {
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;
    typedef std::set<ScaffoldVertex> VertexSet;

    SimpleTransitionGraph graph_;
    ScaffoldVertex source_;
    ScaffoldVertex sink_;
    std::map<VertexSet, size_t> path_cluster_to_weight_;
    vector<ScaffoldVertex> correct_path_;

  SubgraphInfo(const SimpleTransitionGraph &graph,
               const ScaffoldVertex &source,
               const ScaffoldVertex &sink,
               const map<VertexSet, size_t> &path_cluster_to_weight,
               const vector<ScaffoldVertex> &correct_path);

  friend ostream &operator<<(ostream &os, const SubgraphInfo &info);
};

class SubgraphInfoPrinter {
 public:
    void PrintSubgraphInfo(const vector<SubgraphInfo> &info_collection, const string &output_path) const;
};

class PathClusterStatisticsExtractor {
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
    typedef SimpleGraph<scaffold_graph::ScaffoldVertex> SimpleTransitionGraph;

    typedef std::tuple<const Graph&,
                       shared_ptr<cluster_storage::InitialClusterStorage>,
                       shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor>, size_t> PathClusterExtractionParams;

    const conj_graph_pack &gp_;
 public:
    PathClusterStatisticsExtractor(const conj_graph_pack &gp);

    vector<SubgraphInfo> GetAllSubgraphInfo(const ScaffoldGraphStorage &storage);

 private:
    SubgraphInfo GetSubgraphInfo(const SimpleTransitionGraph &graph,
                                 const scaffold_graph::ScaffoldVertex &source,
                                 const scaffold_graph::ScaffoldVertex &sink,
                                 const validation::SimpleTransitionGraphValidator &validator,
                                 const PathClusterExtractionParams &path_cluster_extraction_params) const;

    string RequestCorrectPath(const SimpleTransitionGraph &graph,
                              const scaffold_graph::ScaffoldVertex &source,
                              const scaffold_graph::ScaffoldVertex &sink,
                              const validation::SimpleTransitionGraphValidator &validator) const;

    bool CheckSubgraph(const SimpleTransitionGraph &graph,
                       const scaffold_graph::ScaffoldVertex &source,
                       const scaffold_graph::ScaffoldVertex &sink,
                       const validation::SimpleTransitionGraphValidator &validator) const;

    DECL_LOGGER("PathClusterStatisticsExtractor");
};

}