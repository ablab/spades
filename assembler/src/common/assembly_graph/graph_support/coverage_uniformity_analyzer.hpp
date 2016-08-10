//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#pragma once
#include "assembly_graph/core/graph.hpp"

namespace debruijn_graph {

class CoverageUniformityAnalyzer {
private:
    const Graph& g_;
    const size_t length_bound_;
public:
    CoverageUniformityAnalyzer(const Graph& g, const size_t length_bound): g_(g), length_bound_(length_bound){}
    double CountMedianCoverage() const;
    double UniformityFraction(double allowed_variation, double median_coverage) const;
//first - inside [median* (1 - allowed_variation), median* (1 + allowed_variation)], second-outside
    std::pair<size_t, size_t> TotalLengthsNearMedian(double allowed_variation, double median_coverage) const;
    size_t TotalLongEdgeLength() const;
};
}
