//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "distance_estimation.hpp"
#include "paired_info.hpp"
#include "pair_info_filters.hpp"
#include "smoothing_distance_estimation.hpp"

#include "library/library.hpp"
#include "library/library_data.hpp"

#include "configs/distance_estimation.hpp"

namespace distance_estimation {
    using omnigraph::de::AbstractDistanceEstimator;
    using omnigraph::de::AbstractPairInfoChecker;
    using omnigraph::de::PairedInfoIndexT;
    using omnigraph::de::UnclusteredPairedInfoIndexT;

    class AbstractDistanceEstimatorFabric {
    public:
        typedef AbstractPairInfoChecker<debruijn_graph::Graph> PairInfoChecker;

        virtual std::unique_ptr<omnigraph::de::DistanceEstimator> getDistanceEstimator(const debruijn_graph::Graph &graph,
                                                                      const distance_estimation::UnclusteredPairedInfoIndexT<debruijn_graph::Graph> &index,
                                                                      const omnigraph::de::GraphDistanceFinder &distance_finder,
                                                                      const PairInfoChecker &checker,
                                                                      size_t linkage_distance,
                                                                      size_t max_distance) const = 0;
    };

    class DistanceEstimatorFabric : public AbstractDistanceEstimatorFabric {
    public:
        std::unique_ptr<omnigraph::de::DistanceEstimator> getDistanceEstimator(const debruijn_graph::Graph &graph,
                                                              const distance_estimation::UnclusteredPairedInfoIndexT<debruijn_graph::Graph> &index,
                                                              const omnigraph::de::GraphDistanceFinder &distance_finder,
                                                              const PairInfoChecker &checker,
                                                              size_t linkage_distance,
                                                              size_t max_distance) const override {
            return std::make_unique<omnigraph::de::DistanceEstimator>(graph, index, distance_finder, checker, linkage_distance,
                                                                      max_distance);
        }
    };

    class AbstractScaffoldDistanceEstimatorFabric {
    public:
        typedef AbstractPairInfoChecker<debruijn_graph::Graph> PairInfoChecker;

        virtual std::unique_ptr<omnigraph::de::DistanceEstimator> getDistanceEstimator(const debruijn_graph::Graph &graph,
                                                                      const distance_estimation::UnclusteredPairedInfoIndexT<debruijn_graph::Graph> &histogram,
                                                                      const omnigraph::de::GraphDistanceFinder &dist_finder,
                                                                      const PairInfoChecker &checker,
                                                                      std::function<double(int)> weight_f,
                                                                      size_t linkage_distance, size_t max_distance, size_t threshold,
                                                                      double range_coeff, double delta_coeff,
                                                                      size_t cutoff,
                                                                      size_t min_peak_points,
                                                                      double percentage,
                                                                      double derivative_threshold) const = 0;
    };

    class ScaffoldDistanceEstimatorFabric : public AbstractScaffoldDistanceEstimatorFabric {
    public:
        std::unique_ptr<omnigraph::de::DistanceEstimator> getDistanceEstimator(const debruijn_graph::Graph &graph,
                                                              const distance_estimation::UnclusteredPairedInfoIndexT<debruijn_graph::Graph> &histogram,
                                                              const omnigraph::de::GraphDistanceFinder &dist_finder,
                                                              const PairInfoChecker &checker,
                                                              std::function<double(int)> weight_f,
                                                              size_t linkage_distance, size_t max_distance, size_t threshold,
                                                              double range_coeff, double delta_coeff,
                                                              size_t cutoff,
                                                              size_t min_peak_points,
                                                              double percentage,
                                                              double derivative_threshold) const override {
            return std::unique_ptr<omnigraph::de::DistanceEstimator>(
                    new omnigraph::de::SmoothingDistanceEstimator(graph, histogram, dist_finder, checker, weight_f,
                                                                  linkage_distance, max_distance, threshold,
                                                                  range_coeff, delta_coeff, cutoff, min_peak_points,
                                                                  percentage, derivative_threshold));
        }
    };

    void EstimateWithEstimator(PairedInfoIndexT<debruijn_graph::Graph> &clustered_index,
                               const AbstractDistanceEstimator &estimator,
                               AbstractPairInfoChecker<debruijn_graph::Graph> &checker);

    void RefinePairedInfo(PairedInfoIndexT<debruijn_graph::Graph> &clustered_index,
                          const debruijn_graph::Graph &graph);

    void EstimateScaffoldingDistancesInner(PairedInfoIndexT<debruijn_graph::Graph> &scaffolding_index,
                                           const debruijn_graph::Graph &graph,
                                           const io::SequencingLibrary<debruijn_graph::config::LibraryData> &lib,
                                           const UnclusteredPairedInfoIndexT<debruijn_graph::Graph> &paired_index,
                                           const debruijn_graph::config::smoothing_distance_estimator &ade,
                                           const debruijn_graph::config::distance_estimator &de_config =
                                           debruijn_graph::config::distance_estimator(),
                                           const AbstractScaffoldDistanceEstimatorFabric& distance_estimator_fabric =
                                                   ScaffoldDistanceEstimatorFabric());

    void EstimatePairedDistancesInner(PairedInfoIndexT<debruijn_graph::Graph> &clustered_index,
                                      const debruijn_graph::Graph &graph,
                                      const io::SequencingLibrary<debruijn_graph::config::LibraryData> &lib,
                                      const UnclusteredPairedInfoIndexT<debruijn_graph::Graph> &paired_index,
                                      size_t max_repeat_length = std::numeric_limits<size_t>::max(),
                                      const debruijn_graph::config::distance_estimator &de_config =
                                      debruijn_graph::config::distance_estimator(),
                                      const AbstractDistanceEstimatorFabric& distance_estimator_fabric =
                                              DistanceEstimatorFabric());

    void EstimateScaffoldingDistances(PairedInfoIndexT<debruijn_graph::Graph> &scaffolding_index,
                                           const debruijn_graph::Graph &graph,
                                           const io::SequencingLibrary<debruijn_graph::config::LibraryData> &lib,
                                           const UnclusteredPairedInfoIndexT<debruijn_graph::Graph> &paired_index,
                                           const debruijn_graph::config::smoothing_distance_estimator &ade,
                                           const debruijn_graph::config::distance_estimator &de_config =
                                           debruijn_graph::config::distance_estimator());

    void EstimatePairedDistances(PairedInfoIndexT<debruijn_graph::Graph> &clustered_index,
                                      const debruijn_graph::Graph &graph,
                                      const io::SequencingLibrary<debruijn_graph::config::LibraryData> &lib,
                                      const UnclusteredPairedInfoIndexT<debruijn_graph::Graph> &paired_index,
                                      size_t max_repeat_length = std::numeric_limits<size_t>::max(),
                                      const debruijn_graph::config::distance_estimator &de_config =
                                      debruijn_graph::config::distance_estimator());
}
