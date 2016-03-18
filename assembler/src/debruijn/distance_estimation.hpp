//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {

class DistanceEstimation : public spades::AssemblyStage {
  public:
    DistanceEstimation(bool preliminary = false)
        : AssemblyStage(preliminary ? "Preliminary Distance Estimation" : "Distance Estimation",
                        preliminary ? "distance_estimation_preliminary" : "distance_estimation") {}

    void run(conj_graph_pack &gp, const char*);
};

}

