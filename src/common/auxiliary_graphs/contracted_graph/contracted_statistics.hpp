//***************************************************************************
//* Copyright (c) 2023-2025 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "contracted_graph.hpp"

namespace contracted_graph {
class ContractedStatisticsExtractor {
 public:
    using Graph = debruijn_graph::Graph;

    explicit ContractedStatisticsExtractor(const Graph &assembly_graph);

    size_t CountLoops(const ContractedGraph &graph) const;
    size_t CountNonIsolated(const ContractedGraph &graph) const;
    double GetMeanWeight(const ContractedGraph &graph) const;

    void GetMeanWeights(std::vector<size_t> thresholds, const std::string &output_path) const;

  private:
    const Graph &assembly_graph_;
};
}