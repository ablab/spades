//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "binning.hpp"

namespace bin_stats {

class BinningRefiner {
 public:
    BinningRefiner(const debruijn_graph::Graph& g) : g_(g) {}

    virtual SoftBinsAssignment RefineBinning(const BinStats& bin_stats) const = 0;
 protected:
    const debruijn_graph::Graph& g_;
};
}
