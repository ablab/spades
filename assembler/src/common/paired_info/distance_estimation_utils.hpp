//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "distance_estimation.hpp"
#include "paired_info.hpp"
#include "pair_info_filters.hpp"

#include "library/library.hpp"
#include "library/library_data.hpp"

#include "pipeline/configs/distance_estimation.hpp"

namespace distance_estimation {
using omnigraph::de::AbstractDistanceEstimator;
using omnigraph::de::AbstractPairInfoChecker;
using omnigraph::de::PairedInfoIndexT;
using omnigraph::de::UnclusteredPairedInfoIndexT;

void EstimateWithEstimator(PairedInfoIndexT<debruijn_graph::Graph> &clustered_index,
                           const AbstractDistanceEstimator &estimator,
                           AbstractPairInfoChecker<debruijn_graph::Graph> &checker);

void RefinePairedInfo(PairedInfoIndexT<debruijn_graph::Graph>& clustered_index,
                      const debruijn_graph::Graph& graph);

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
