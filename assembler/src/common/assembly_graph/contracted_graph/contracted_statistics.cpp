#include "contracted_statistics.hpp"
#include "contracted_graph_builder.hpp"
namespace contracted_graph {
size_t ContractedStatisticsExtractor::CountLoops(const ContractedGraph &graph) const {
    size_t result = 0;
    for (const auto &vertex: graph) {
        for (auto it = graph.out_begin(vertex); it != graph.out_end(vertex); ++it) {
            auto next_vertex = it->first;
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
        if (graph.getOutDegree(vertex) > 0 or graph.getInDegree(vertex) > 0) {
            total_capacity += graph.capacity(vertex);
            ++non_isolated;
        }
    }
    INFO("Total capacity: " << total_capacity << ", Non isolated: " << non_isolated);
    return static_cast<double>(total_capacity) / static_cast<double>(non_isolated);
}
void ContractedStatisticsExtractor::GetMeanWeights(vector<size_t> thresholds, const string &output_path) const {
    std::map<size_t, double> threshold_to_mean_weight;
    for (const size_t threshold: thresholds) {
        INFO("Constructing graph for " << threshold);
        auto length_predicate = [this, threshold](const EdgeId &edge) {
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
    ofstream fout(output_path);
    for (const auto &entry: threshold_to_mean_weight) {
        fout << entry.first << " " << entry.second << "\n";
    }
}
ContractedStatisticsExtractor::ContractedStatisticsExtractor(const Graph &assembly_graph) : assembly_graph_(
    assembly_graph) {}
size_t ContractedStatisticsExtractor::CountNonIsolated(const ContractedGraph &graph) const {
    size_t result = 0;
    for (const auto &vertex: graph) {
        if (graph.getInDegree(vertex) > 0 or graph.getOutDegree(vertex) > 0) {
            ++result;
        }
    }
    return result;
}
}
