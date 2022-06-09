#include "common/assembly_graph/core/graph.hpp"
#include "common/assembly_graph/core/basic_graph_stats.hpp"
#include "common/assembly_graph/stats/statistics.hpp"
#include "common/toolchain/subgraph_utils.hpp"
#include "common/assembly_graph/components/graph_component.hpp"
#include "common/assembly_graph/components/connected_component.hpp"

#include <iostream>

/*
void CountDegree(const debruijn_graph::Graph& g) {
    size_t max_degree = 0;
    size_t min_degree = 1e9;
    double average_degree = 0;

    for (auto it = g.SmartVertexBegin(); !it.IsEnd(); ++it) {

    }

    std::cout << "Average degree of the sequences : ";
    std::cout << "Min degree of the sequences : ";
    std::cout << "Max degree of the sequences : ";
}
*/

void PrintStat(const debruijn_graph::Graph& g, const std::vector<double>& n_percentiles,
               const std::vector<double>& median_length_percentiles,
               const std::vector<double>& sequences_coverage_percentiles) {
    std::cout << "k : " << g.k() << "\n";

    size_t n_seq = g.e_size();
    std::cout << "Number of sequences : " << n_seq << "\n";

    size_t n_links = g.size();
    std::cout << "Number of links : " << n_links << "\n";

    omnigraph::CumulativeLengthCounter cumulative_length_counter(g);
    std::cout << "Total sequences length : " << cumulative_length_counter.Count() << "\n";

    std::vector<size_t> sequences_lengths;
    double average_length = 0;
    for (auto it = g.ConstEdgeBegin(); !it.IsEnd(); ++it) {
        size_t sequence_length = g.EdgeNucls(*it).str().length();
        sequences_lengths.push_back(sequence_length);
        average_length += sequence_length;
    }
    std::cout << "Average sequences length : " << average_length / n_seq << "\n";
    std::sort(sequences_lengths.begin(), sequences_lengths.end());
    std::cout << "Shortest sequence : " << sequences_lengths[0] << "\n";
    for (double perc_ : median_length_percentiles) {
        size_t perc_end_idx = sequences_lengths.size() * perc_ / 100;
        size_t median_idx = perc_end_idx / 2;
        std::cout << "Median length for " << perc_ << " percentile of sequences : " <<
                                                        sequences_lengths[median_idx] << "\n";
    }
    std::cout << "Longest sequence : " << *sequences_lengths.rbegin() << "\n";

    omnigraph::GraphComponent<debruijn_graph::Graph> graph_component = omnigraph::GraphComponent<debruijn_graph::Graph>::WholeGraph(g);
    std::vector<debruijn_graph::EdgeId> dead_ends = toolchain::CollectDeadEnds(graph_component);
    size_t n_dead_ends = dead_ends.size();
    std::cout << "Number of sequences with dead ends : " << n_dead_ends << "\n";
    std::cout << "Percentage of sequences with dead ends : " << n_dead_ends * 100.0 / n_seq ;

    debruijn_graph::ConnectedComponentCounter connected_component_counter(g);
    connected_component_counter.CalculateComponents();
    size_t n_components = 0;
    for (const auto& edge : connected_component_counter.component_ids_) {
        n_components = std::max(n_components, edge.second);
    }

    std::cout << "Number of connected components : " << n_components << "\n";

    size_t largest_component_id = 0;
    size_t largest_component_size = 0;
    for (const auto& component : connected_component_counter.component_total_len_) {
        if (component.second > largest_component_size) {
            largest_component_id = component.first;
            largest_component_size = component.second;
        }
    }

    std::cout << "The total number of base pairs in the largest connected component : " << largest_component_size << "\n";

    for (double perc_ : n_percentiles) {
        std::cout << "N" << perc_ << " : " << omnigraph::Nx/*<debruijn_graph::Graph>*/(g, perc_) << "\n";
    }


    std::vector<double> average_sequences_coverages;
    for (auto it = g.ConstEdgeBegin(); !it.IsEnd(); ++it) {
        average_sequences_coverages.push_back(g.coverage(*it));
    }
    std::sort(average_sequences_coverages.begin(), average_sequences_coverages.end());
    for (double perc_ : sequences_coverage_percentiles) {
        size_t idx = average_sequences_coverages.size() * perc_ / 100;
        std::cout << perc_ <<" percentile sequences coverage : " << average_sequences_coverages[idx] << "\n";
    }

    //std::cout << "Estimated genome length : ";
}