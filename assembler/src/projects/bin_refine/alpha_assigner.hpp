//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "binning.hpp"
#include "id_map.hpp"

#include "assembly_graph/core/graph.hpp"

namespace bin_stats {

using AlphaAssignment = adt::id_map<double, debruijn_graph::EdgeId>;

class AlphaAssigner {
  public:
    explicit AlphaAssigner(const debruijn_graph::Graph &g)
        : g_(g) {}
    virtual ~AlphaAssigner() = default;

    virtual AlphaAssignment GetAlphaAssignment(const Binning &bin_stats) const = 0;

  protected:
    //fixme duplication with LabelsPropagation
    SoftBinsAssignment InitLabels(const Binning& bin_stats) const;

    const debruijn_graph::Graph &g_;
};

class PropagationAssigner : public AlphaAssigner {
  public:
    explicit PropagationAssigner(const debruijn_graph::Graph &g)
        : AlphaAssigner(g) {}

    AlphaAssignment GetAlphaAssignment(const Binning &bin_stats) const override;

  private:
    using AlphaAssigner::g_;
};

}