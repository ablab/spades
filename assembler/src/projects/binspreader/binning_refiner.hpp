//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "binning.hpp"

namespace binning {
class LinkIndex;
}

namespace bin_stats {

class BinningRefiner {
 public:
    explicit BinningRefiner(const debruijn_graph::Graph& g,
                            const binning::LinkIndex &links)
            : g_(g), links_(links) {}
    virtual ~BinningRefiner() = default;

    virtual SoftBinsAssignment RefineBinning(const SoftBinsAssignment &state) const = 0;
 protected:
    const debruijn_graph::Graph& g_;
    const binning::LinkIndex &links_;
};
}
