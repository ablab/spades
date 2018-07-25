#pragma once

#include <ostream>
#include "common/modules/path_extend/read_cloud_path_extend/validation/transition_subgraph_validation.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/path_cluster_helper.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/simple_graph.hpp"
#include "common/barcode_index/cluster_storage_extractor.hpp"

namespace path_extend {

struct SubgraphInfo {
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;
    typedef std::set<ScaffoldVertex> VertexSet;
    typedef std::map<ScaffoldVertex, std::map<ScaffoldVertex, size_t>> ScaffoldEdgeMap;


    SimpleTransitionGraph graph_;
    ScaffoldVertex source_;
    ScaffoldVertex sink_;
    PathClusterStorage path_cluster_to_weight_;
    vector<VertexSet> final_clusters_;
    vector<vector<ScaffoldVertex>> resulting_paths_;
    vector<vector<ScaffoldVertex>> all_paths_;
    vector<ScaffoldVertex> correct_path_;
    std::map<SubgraphInfo::ScaffoldVertex, string> id_map_;
    std::map<ScaffoldVertex, double> vertex_to_cov_;
    std::map<ScaffoldVertex, size_t> vertex_to_len_;
    ScaffoldEdgeMap scaffold_edge_to_dist_;


  SubgraphInfo(const SimpleTransitionGraph &graph,
               const ScaffoldVertex &source,
               const ScaffoldVertex &sink,
               const PathClusterStorage &path_cluster_to_weight,
               const vector<VertexSet> &final_clusters,
               const vector<vector<ScaffoldVertex>> &resulting_paths,
               const vector<vector<ScaffoldVertex>> &all_paths,
               const vector<ScaffoldVertex> &correct_path,
               const std::map<SubgraphInfo::ScaffoldVertex, string> &id_map,
               const std::map<ScaffoldVertex, double> &vertex_to_cov,
               const std::map<ScaffoldVertex, size_t> &vertex_to_len,
               const ScaffoldEdgeMap &scaffold_edge_to_dist);

  friend ostream &operator<<(ostream &os, const SubgraphInfo &info);
};

class SubgraphInfoPrinter {
 public:
    void PrintSubgraphInfo(const vector<SubgraphInfo> &info_collection, const string &output_path) const;
};

struct PathClusterExtractionParams {
  const Graph &g_;
  shared_ptr<cluster_storage::InitialClusterStorage> init_cluster_storage_;
  shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
  size_t linkage_distance_;

  PathClusterExtractionParams(const Graph &g,
                              shared_ptr<cluster_storage::InitialClusterStorage> init_cluster_storage,
                              shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                              size_t linkage_distance);
};

class PathClusterStatisticsExtractor {
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;
    typedef std::map<ScaffoldVertex, string> IdMap;
    typedef std::map<ScaffoldVertex, std::map<ScaffoldVertex, size_t>> ScaffoldEdgeMap;

    const conj_graph_pack &gp_;
 public:
    PathClusterStatisticsExtractor(const conj_graph_pack &gp);

    vector<SubgraphInfo> GetAllSubgraphInfo(const ScaffoldGraphStorage &storage);

 private:
    SubgraphInfo GetSubgraphInfo(const SimpleTransitionGraph &graph,
                                 const ScaffoldVertex &source,
                                 const ScaffoldVertex &sink,
                                 const ScaffoldEdgeMap &scaffold_edge_to_len,
                                 const validation::SimpleTransitionGraphValidator &validator,
                                 const PathClusterExtractionParams &path_cluster_extraction_params,
                                 bool reference_validation_on) const;

    string RequestCorrectPath(const SimpleTransitionGraph &graph,
                              const ScaffoldVertex &source,
                              const ScaffoldVertex &sink,
                              const validation::SimpleTransitionGraphValidator &validator) const;

    bool CheckSubgraph(const SimpleTransitionGraph &graph,
                       const ScaffoldVertex &source,
                       const ScaffoldVertex &sink,
                       const validation::SimpleTransitionGraphValidator &validator,
                       bool reference_validation_on) const;

    IdMap GetIdMap(const SimpleTransitionGraph &graph, const vector<ScaffoldVertex> &correct_path_result) const;

    ScaffoldEdgeMap ConstructLengthMap(const SimpleTransitionGraph &transition_graph, const ScaffoldGraph &graph) const;

    DECL_LOGGER("PathClusterStatisticsExtractor");
};

}