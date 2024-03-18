//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "alpha_assigner.hpp"
#include "binning.hpp"
#include "link_index.hpp"

#include <filesystem>

namespace bin_stats {

class AlphaPropagator {
  public:
    using EdgeId = debruijn_graph::EdgeId;
    using Graph = debruijn_graph::Graph;
    using BinId = Binning::BinId;

    static constexpr BinId UNBINNED = BinId(0);
    static constexpr BinId BINNED = BinId(1);

    AlphaPropagator(const debruijn_graph::Graph &g,
                    const binning::LinkIndex &links,
                    double metaalpha,
                    double eps, unsigned niter,
                    size_t length_threshold,
                    size_t distance_threshold,
                    const std::filesystem::path &debug_path)
        : g_(g),
          links_(links),
          metaalpha_(metaalpha),
          eps_(eps),
          niter_(niter),
          length_threshold_(length_threshold),
          distance_bound_(distance_threshold),
          debug_path_(debug_path) {}

    AlphaAssignment GetAlphaMask(const Binning &bin_stats) const;
  private:
    SoftBinsAssignment ConstructBinningMask(const SoftBinsAssignment &origin_state) const;

    const Graph &g_;
    const binning::LinkIndex &links_;
    double metaalpha_;
    double eps_;
    unsigned niter_;
    size_t length_threshold_;
    size_t distance_bound_;
    std::filesystem::path debug_path_;
};

}
