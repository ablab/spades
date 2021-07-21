//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "alpha_assigner.hpp"
#include "binning.hpp"
#include "link_index.hpp"

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
                    double eps, 
                    const std::string &debug_path)
        : g_(g), links_(links), metaalpha_(metaalpha), eps_(eps), debug_path_(debug_path) {}

    AlphaAssignment GetAlphaMask(const Binning &bin_stats) const;
  private:
    SoftBinsAssignment ConstructBinningMask(const SoftBinsAssignment &origin_state) const;

    const Graph &g_;
    const binning::LinkIndex links_;
    double metaalpha_;
    double eps_;
    std::string debug_path_;
};

}