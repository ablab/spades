//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/path_extend/path_extender.hpp"

namespace path_extend {
namespace read_cloud {

struct SearchParams {
  SearchParams();
  SearchParams(size_t max_path_growing_iterations,
               size_t max_paths_to_process,
               size_t max_edge_visits);
  const size_t max_path_growing_iterations;
  const size_t max_paths_to_process;
  const size_t max_edge_visits;
};

struct SearchingExtenderParams {
  explicit SearchingExtenderParams(const ScaffoldingUniqueEdgeStorage &unique_storage);
  SearchingExtenderParams(const ScaffoldingUniqueEdgeStorage &unique_storage,
                          size_t insert_size,
                          bool investigate_short_loops,
                          bool use_short_loop_cov_resolver,
                          double weight_threshold,
                          size_t length_bound);
  const ScaffoldingUniqueEdgeStorage &unique_storage;
  size_t insert_size;
  bool investigate_short_loops;
  bool use_short_loop_cov_resolver;
  double weight_threshold;
  size_t length_bound;
};

struct ChooserConstructionParams {
  explicit ChooserConstructionParams(const config::dataset::Library &lib);
  ChooserConstructionParams(const config::dataset::Library &lib,
                            size_t lib_index,
                            bool is_coverage_aware,
                            double lib_cov,
                            double single_threshold,
                            bool normalize_weight,
                            double weight_threshold,
                            double priority_coeff);

  const debruijn_graph::config::dataset::Library &lib;
  size_t lib_index;
  bool is_coverage_aware;
  double lib_cov;
  double single_threshold;
  bool normalize_weight;
  double weight_threshold;
  double priority_coeff;
};

struct ReadCloudSearchParameterPack {
  ChooserConstructionParams chooser_params;
  SearchParams search_params;
  SearchingExtenderParams searching_extender_params;
};


class DefaultExtenderParamsConstructor {
  public:
    using GraphPack = graph_pack::GraphPack;

    DefaultExtenderParamsConstructor(const graph_pack::GraphPack &gp,
                                     const config::dataset &dataset_info,
                                     const ScaffoldingUniqueEdgeStorage &unique_storage);

    SearchingExtenderParams ConstructExtenderParams(size_t lib_index, double weight_threshold) const;

  private:
    const graph_pack::GraphPack &gp_;
    const config::dataset &dataset_info_;
    const ScaffoldingUniqueEdgeStorage &unique_storage_;
};
//
//class ExtenderSearcher {
//  public:
//
//    ExtenderSearcher(const conj_graph_pack &gp, std::shared_ptr<ExtensionChooser> extension_chooser,
//                     const SearchParams &params,
//                     const SearchingExtenderParams &extender_params, size_t length_bound);
//
//    bool IsReachable(VertexId target_vertex, BidirectionalPath *start) const;
//
//  private:
//    bool SearchTargetUsingExtenders(std::shared_ptr<path_extend::QueueContainer> paths_container,
//                                    GraphCoverageMap &cover_map,
//                                    const VertexId &target_vertex) const;
//
//    std::shared_ptr<SearchingMultiExtender> ConstructSearchingExtender(std::shared_ptr<ExtensionChooser> extension_chooser,
//                                                                       std::shared_ptr<QueueContainer> paths_container,
//                                                                       GraphCoverageMap &cover_map,
//                                                                       UsedUniqueStorage &used_unique_storage) const;
//
//    const conj_graph_pack &gp_;
//    std::shared_ptr<ExtensionChooser> extension_chooser_;
//    SearchParams search_params_;
//    SearchingExtenderParams extender_params_;
//    size_t length_bound_;
//
//    DECL_LOGGER("ExtenderSearcher");
//};
}
}