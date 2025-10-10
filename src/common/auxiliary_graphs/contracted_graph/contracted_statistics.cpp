//***************************************************************************
//* Copyright (c) 2023-2025 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "contracted_statistics.hpp"

#include "contracted_graph_builder.hpp"

namespace contracted_graph {
size_t ContractedStatisticsExtractor::CountLoops(const ContractedGraph &graph) const {
    size_t result = 0;
    for (const auto &vertex: graph) {
        for (const auto &entry: graph.OutcomingEntries(vertex)) {
            auto next_vertex = entry.first;
            if (vertex == next_vertex) {
                ++result;
            }
        }
    }
    return result;
}
double ContractedStatisticsExtractor::GetMeanWeight(const ContractedGraph &graph) const {
    size_t total_capacity = 0;
    size_t non_isolated = 0;
    for (const auto &vertex: graph) {
        if (graph.GetOutDegree(vertex) > 0 or graph.GetInDegree(vertex) > 0) {
            total_capacity += graph.GetCapacity(vertex);
            ++non_isolated;
        }
    }
    INFO("Total capacity: " << total_capacity << ", Non isolated: " << non_isolated);
    return static_cast<double>(total_capacity) / static_cast<double>(non_isolated);
}
void ContractedStatisticsExtractor::GetMeanWeights(std::vector<size_t> thresholds,
                                                   const std::string &output_path) const {
    std::map<size_t, double> threshold_to_mean_weight;
    for (const size_t threshold: thresholds) {
        INFO("Constructing graph for " << threshold);
        auto length_predicate = [this, threshold](debruijn_graph::EdgeId edge) {
          return assembly_graph_.length(edge) >= threshold;
        };
        DBGContractedGraphFactory factory(assembly_graph_, length_predicate);
        factory.Construct();
        INFO("Constructed graph");
        auto graph = factory.GetGraph();
        double mean_weight = GetMeanWeight(*graph);
        INFO("Mean weight " << mean_weight);
        threshold_to_mean_weight.insert({threshold, mean_weight});
    }
    std::ofstream fout(output_path);
    for (const auto &entry: threshold_to_mean_weight) {
        fout << entry.first << " " << entry.second << "\n";
    }
}
ContractedStatisticsExtractor::ContractedStatisticsExtractor(const Graph &assembly_graph) :
    assembly_graph_(assembly_graph) {}
size_t ContractedStatisticsExtractor::CountNonIsolated(const ContractedGraph &graph) const {
    size_t result = 0;
    for (const auto &vertex: graph) {
        if (graph.GetInDegree(vertex) > 0 or graph.GetOutDegree(vertex) > 0) {
            ++result;
        }
    }
    return result;
}
}
