//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/paired_info/distance_estimation_utils.hpp"
#include "distance_estimation.hpp"

namespace distance_estimation {
    using omnigraph::de::DistanceEstimator;
    using omnigraph::de::DistanceEstimatorMPI;

    class MPIDistanceEstimatorFabric : public AbstractDistanceEstimatorFabric {
    public:
        std::unique_ptr<DistanceEstimator> getDistanceEstimator(const debruijn_graph::Graph &graph,
                                                                const distance_estimation::UnclusteredPairedInfoIndexT<debruijn_graph::Graph> &index,
                                                                const omnigraph::de::GraphDistanceFinder &distance_finder,
                                                                const PairInfoChecker &checker,
                                                                size_t linkage_distance,
                                                                size_t max_distance) const override {
            auto estimator_base = std::make_unique<DistanceEstimator>(graph, index, distance_finder, checker,
                                                                      linkage_distance, max_distance);
            return std::unique_ptr<DistanceEstimator>(new DistanceEstimatorMPI(graph, index,
                                                                               distance_finder,
                                                                               checker,
                                                                               linkage_distance,
                                                                               max_distance,
                                                                               std::move(estimator_base)));
        }
    };

    class MPIScaffoldDistanceEstimatorFabric : public AbstractScaffoldDistanceEstimatorFabric {
    public:
        std::unique_ptr<omnigraph::de::DistanceEstimator> getDistanceEstimator(const debruijn_graph::Graph &graph,
                                                                               const distance_estimation::UnclusteredPairedInfoIndexT<debruijn_graph::Graph> &histogram,
                                                                               const omnigraph::de::GraphDistanceFinder &dist_finder,
                                                                               const PairInfoChecker &checker,
                                                                               std::function<double(int)> weight_f,
                                                                               size_t linkage_distance,
                                                                               size_t max_distance, size_t threshold,
                                                                               double range_coeff, double delta_coeff,
                                                                               size_t cutoff,
                                                                               size_t min_peak_points,
                                                                               double percentage,
                                                                               double derivative_threshold) const override {
            auto estimator_base = std::unique_ptr<omnigraph::de::DistanceEstimator>(
                new omnigraph::de::SmoothingDistanceEstimator(graph, histogram, dist_finder, checker, weight_f,
                                                                  linkage_distance, max_distance, threshold,
                                                                  range_coeff, delta_coeff, cutoff, min_peak_points,
                                                                  percentage, derivative_threshold));

            return std::unique_ptr<DistanceEstimator>(new DistanceEstimatorMPI(graph, histogram,
                                                                               dist_finder,
                                                                               checker,
                                                                               linkage_distance,
                                                                               max_distance,
                                                                               std::move(estimator_base)));
        }
    };

    void EstimateScaffoldingDistancesMPI(PairedInfoIndexT<debruijn_graph::Graph> &scaffolding_index,
                                         const debruijn_graph::Graph &graph,
                                         const io::SequencingLibrary<debruijn_graph::config::LibraryData> &lib,
                                         const UnclusteredPairedInfoIndexT<debruijn_graph::Graph> &paired_index,
                                         const debruijn_graph::config::smoothing_distance_estimator &ade,
                                         const debruijn_graph::config::distance_estimator &de_config =
                                         debruijn_graph::config::distance_estimator());

    void EstimatePairedDistancesMPI(PairedInfoIndexT<debruijn_graph::Graph> &clustered_index,
                                    const debruijn_graph::Graph &graph,
                                    const io::SequencingLibrary<debruijn_graph::config::LibraryData> &lib,
                                    const UnclusteredPairedInfoIndexT<debruijn_graph::Graph> &paired_index,
                                    size_t max_repeat_length = std::numeric_limits<size_t>::max(),
                                    const debruijn_graph::config::distance_estimator &de_config =
                                    debruijn_graph::config::distance_estimator());
}
