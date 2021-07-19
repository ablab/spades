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

class AlphaPropagator : public AlphaAssigner {
  public:
    using EdgeId = debruijn_graph::EdgeId;
    using Graph = debruijn_graph::Graph;

    AlphaPropagator(const debruijn_graph::Graph &g,
                    const binning::LinkIndex &links,
                    double metaalpha)
        : AlphaAssigner(g), links_(links), metaalpha_(metaalpha) {}

    AlphaAssignment GetAlphaAssignment(const Binning &bin_stats) const override;
  private:
    SoftBinsAssignment ConstructBinningMask(const SoftBinsAssignment &origin_state) const;

    const binning::LinkIndex links_;
    double metaalpha_;
};

}