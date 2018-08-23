#pragma once
#include "contracted_graph.hpp"


namespace contracted_graph {
class ContractedStatisticsExtractor {
    const Graph& assembly_graph_;
 public:
    explicit ContractedStatisticsExtractor(const Graph &assembly_graph);

    size_t CountLoops(const ContractedGraph &graph) const;
    size_t CountNonIsolated(const ContractedGraph &graph) const;
    double GetMeanWeight(const ContractedGraph &graph) const;

    void GetMeanWeights(vector<size_t> thresholds, const string &output_path) const;
};
}