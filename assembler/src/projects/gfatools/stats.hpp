#include "common/assembly_graph/core/graph.hpp"
#include "common/assembly_graph/core/basic_graph_stats.hpp"
#include "common/assembly_graph/stats/statistics.hpp"
#include "common/toolchain/subgraph_utils.hpp"
#include "common/assembly_graph/components/graph_component.hpp"
#include "common/assembly_graph/components/connected_component.hpp"
#include "io/binary/binary.hpp"

#include "llvm/Support/YAMLTraits.h"
#include "llvm/Support/FileSystem.h"

#include <iostream>
#include <map>

struct StatInfo {

    struct GraphSize {
        size_t k;
        size_t n_seq;
        size_t n_links;
        uint64_t total_length;
        uint64_t total_length_no_overlaps;
    };


    struct SequencesLengthInfo {
        std::map<double, size_t> nx;
        double average_length;
        size_t longest_sequence;
        size_t shortest_sequence;
        std::map<double, double> percentile_length;
    };

    struct GraphConnectivityInfo {
        size_t n_dead_ends;
        double dead_ends_percentage;
        size_t n_connected_components;
        uint64_t largest_component_length;
        size_t max_degree;
        double average_degree;
    };

    struct CoverageInfo {
        std::map<double, double> percentile_coverage;
        uint64_t estimated_length;
    };
    GraphSize graph_size;
    SequencesLengthInfo length_info;
    GraphConnectivityInfo connectivity_info;
    CoverageInfo coverage_info;
};

namespace llvm { namespace yaml {

template<> struct MappingTraits<StatInfo> {
    static void mapping(IO& io, StatInfo& info) {
        io.mapRequired("Graph size", info.graph_size);
        io.mapRequired("Sequences length info", info.length_info);
        io.mapRequired("Graph connectivity info", info.connectivity_info);
        io.mapRequired("Coverage info", info.coverage_info);
    }
};
 template<> struct MappingTraits<StatInfo::GraphSize> {
     static void mapping(IO& io, StatInfo::GraphSize& info) {
         io.mapRequired("Sequences overlap", info.k);
         io.mapRequired("Number of sequences", info.n_seq);
         io.mapRequired("Number of links", info.n_links);
         io.mapRequired("Total length", info.total_length);
         io.mapRequired("Total length no overlaps", info.total_length_no_overlaps);
     }
 };

 template<> struct MappingTraits<StatInfo::SequencesLengthInfo> {
     static void mapping(IO& io, StatInfo::SequencesLengthInfo& info) {
         io.mapRequired("N stat for percentiles", info.nx);
         io.mapRequired("Average sequence length", info.average_length);
         io.mapRequired("Longest sequence", info.longest_sequence);
         io.mapRequired("Shortest sequence", info.shortest_sequence);
         io.mapRequired("Length for percentiles", info.percentile_length);
     }
 };

 template<> struct MappingTraits<StatInfo::GraphConnectivityInfo> {
     static void mapping(IO& io, StatInfo::GraphConnectivityInfo& info) {
         io.mapRequired("Number of dead ends", info.n_dead_ends);
         io.mapRequired("Percentage of dead ends", info.dead_ends_percentage);
         io.mapRequired("Number of connected components", info.n_connected_components);
         io.mapRequired("Largest component", info.largest_component_length);
         io.mapRequired("Max degree", info.max_degree);
         io.mapRequired("Average degree", info.average_degree);
     }
 };

 template<> struct MappingTraits<StatInfo::CoverageInfo> {
     static void mapping(IO& io, StatInfo::CoverageInfo& info) {
         io.mapRequired("Coverage for percentile", info.percentile_coverage);
         io.mapRequired("Estimated genome length", info.estimated_length);
     }
 };

} }

template<class T>
void WriteToYaml(T &t, const std::filesystem::path& filename) {
    std::ofstream ofs(filename, std::ios::binary);
    WriteToYaml(t, ofs);
}

template<class T>
void WriteToYaml(T &t, std::ostream& os) {
    std::string buffer;
    llvm::raw_string_ostream raw_os(buffer);
    llvm::yaml::Output yout(raw_os);
    yout << t;
    raw_os.str();
    io::binary::BinWrite(os, buffer);
}

void GetSequencesNumber(const debruijn_graph::Graph& g, size_t& n_seq, StatInfo& info) {
    for (debruijn_graph::EdgeId e: g.canonical_edges()) {
        n_seq++;
    }
    info.graph_size.n_seq = n_seq;
}


void GetLinksNumber(const debruijn_graph::Graph& g, StatInfo& info) {
    size_t n_links = 0;
    for (debruijn_graph::VertexId v : g.canonical_vertices()) {
        size_t in =  g.IncomingEdgeCount(v);
        size_t out = g.OutgoingEdgeCount(v);
        n_links += in * out;
    }
    info.graph_size.n_links = n_links;
}


void GetDeadEndsInfo(const debruijn_graph::Graph& g, StatInfo& info) {
    omnigraph::GraphComponent<debruijn_graph::Graph> graph_component = omnigraph::GraphComponent<debruijn_graph::Graph>::WholeGraph(g);
    std::vector<debruijn_graph::EdgeId> dead_ends = toolchain::CollectDeadEnds(graph_component);
    size_t n_dead_ends = dead_ends.size();
    info.connectivity_info.n_dead_ends = n_dead_ends;
    info.connectivity_info.dead_ends_percentage = n_dead_ends * 100.0 / g.e_size();
}


template<typename T>
double GetPercentile(const std::vector<T>& v, double perc) {
    if (perc < 0 || perc > 100) {
        return -1;
    }
    if (perc == 0) {
        return v[0];
    }
    if (perc == 100) {
        return v[v.size() - 1];
    }
    double perc_idx = v.size() * perc / 100;
    size_t integer_part = floor(perc_idx);
    double fractional_part = perc_idx - integer_part;
    if (fractional_part > 1e-9) {
        return v[integer_part];
    }
    else {
        return 1.0 * (v[integer_part - 1] + v[integer_part]) / 2.0;
    }
}


void GetSequencesLengthInfo(const debruijn_graph::Graph& g,
                            const std::vector<double>& percentiles, size_t n_seq, StatInfo& info) {
    std::vector<size_t> sequences_lengths;
    uint64_t total_length = 0;
    uint64_t total_length_no_overlaps = 0;
    for (debruijn_graph::EdgeId e: g.edges()) {
        size_t sequence_length = g.EdgeNucls(e).size();
        sequences_lengths.push_back(sequence_length);
        total_length += sequence_length;
    }

    uint64_t doubles = 0;
    for (debruijn_graph::VertexId v : g.vertices()) {
        doubles += (g.IncomingEdgeCount(v) + g.OutgoingEdgeCount(v) - 1);
    }
    for (debruijn_graph::EdgeId e : g.edges()) {
        if (g.EdgeStart(e) == g.EdgeEnd(e)) {
            doubles--;
        }
    }
    info.graph_size.total_length = total_length / 2;
    info.graph_size.total_length_no_overlaps = (total_length - g.k() * doubles) / 2;

    info.length_info.average_length = 1.0 * total_length / (2 * n_seq);

    std::sort(sequences_lengths.begin(), sequences_lengths.end());
    info.length_info.shortest_sequence = sequences_lengths[0];
    for (double perc_ : percentiles) {
        double res = GetPercentile(sequences_lengths, perc_);
        if (res == -1) {
            continue;
        }
        info.length_info.percentile_length[perc_] = res;
    }

    info.length_info.longest_sequence = *sequences_lengths.rbegin();
}


void GetConnectedComponentsInfo(const debruijn_graph::Graph& g, StatInfo& info) {
    debruijn_graph::ConnectedComponentCounter connected_component_counter(g);
    connected_component_counter.CalculateComponents();

    size_t n_components = 0;
    for (const auto& edge : connected_component_counter.component_ids_) {
        n_components = std::max(n_components, edge.second + 1);
    }
    info.connectivity_info.n_connected_components = n_components;

    std::map<size_t, size_t> components_length;
    for (debruijn_graph::EdgeId e : g.canonical_edges()) {
        components_length[connected_component_counter.component_ids_[e]] += g.EdgeNucls(e).size();
    }

    size_t largest_component_size = 0;
    for (const auto& component : components_length) {
        largest_component_size = std::max(largest_component_size, component.second);
    }
    info.connectivity_info.largest_component_length = largest_component_size;
}


size_t Nx(const debruijn_graph::Graph& g, double percent) {
    size_t sum_edge_length = 0;
    std::vector<size_t> lengths;
    for (debruijn_graph::EdgeId e : g.canonical_edges()) {
        lengths.push_back(g.EdgeNucls(e).size());
        sum_edge_length += g.EdgeNucls(e).size();
    }
    std::sort(lengths.rbegin(), lengths.rend());
    double len_perc = percent * 0.01 * sum_edge_length;
    for (size_t i = 0; i < lengths.size(); i++) {
        if (double(lengths[i]) >= len_perc)
            return lengths[i];
        else
            len_perc -= (double) lengths[i];
    }
    return 0;
}


void GetPercentileCoverage(const debruijn_graph::Graph& g, const std::vector<double>& percentiles,
                           double& median_coverage, StatInfo& info) {
    std::vector<std::pair<double, size_t>> edges;
    uint64_t total_length = 0;
    for (debruijn_graph::EdgeId e: g.canonical_edges()) {
        edges.push_back(std::make_pair(g.coverage(e), g.EdgeNucls(e).size()));
        total_length += g.EdgeNucls(e).size();
    }
    std::sort(edges.begin(), edges.end());
    uint64_t cumulative_length = 0;
    size_t perc_ind = 0;
    bool median_is_counted = false;
    for (auto e : edges) {
        cumulative_length += e.second;
        while (perc_ind < percentiles.size() && cumulative_length >= total_length * percentiles[perc_ind] / 100.0) {
            info.coverage_info.percentile_coverage[perc_ind++] = e.first;
        }
        if (!median_is_counted && cumulative_length >= total_length * 0.5) {
            median_coverage = e.first;
            median_is_counted = true;
        }
    }
}


void GetEstimatedGenomeLength(const debruijn_graph::Graph& g, double median_coverage, StatInfo& info) {
    uint64_t estimated_length = 0;
    std::map<debruijn_graph::VertexId, size_t> vertexes_degrees;
    for (debruijn_graph::VertexId v : g.vertices()) {
        vertexes_degrees[v] = g.IncomingEdgeCount(v) + g.OutgoingEdgeCount(v);
    }
    for (debruijn_graph::EdgeId e: g.edges()) {
        size_t edge_length = g.EdgeNucls(e).size();
        if (edge_length > g.k() && vertexes_degrees[g.EdgeStart(e)] > 1) {
            edge_length -= g.k();
            vertexes_degrees[g.EdgeStart(e)]--;
        }
        if (g.EdgeStart(e) != g.EdgeEnd(e) && edge_length > g.k() && vertexes_degrees[g.EdgeEnd(e)] > 1) {
            edge_length -= g.k();
            vertexes_degrees[g.EdgeEnd(e)]--;
        }
        double coverage = g.coverage(e);
        size_t relative_coverage = floor(coverage / median_coverage);
        estimated_length += edge_length * relative_coverage;
    }
    info.coverage_info.estimated_length = estimated_length / 2;
}


void GetSequencesDegreeInfo(const debruijn_graph::Graph& g, size_t n_seq, StatInfo& info) {
    size_t max_degree = 0;
    double average_degree = 0;
    for (debruijn_graph::EdgeId e : g.canonical_edges()) {
        debruijn_graph::VertexId start_v = g.EdgeStart(e);
        debruijn_graph::VertexId end_v = g.EdgeEnd(e);
        max_degree = std::max({max_degree, g.IncomingEdgeCount(start_v), g.OutgoingEdgeCount(end_v)});
        average_degree += g.IncomingEdgeCount(start_v);
        average_degree +=  g.OutgoingEdgeCount(end_v);
    }

    info.connectivity_info.max_degree = max_degree;
    info.connectivity_info.average_degree = average_degree / (2 * n_seq);
}

void PrintStat(const StatInfo& info) {
    std::cout.precision(3);
    std::cout << std::fixed;

    std::cout << "Graph size\n";
    std::cout << "Sequences overlap: " << info.graph_size.k << "\n";
    std::cout << "Number of sequences: " << info.graph_size.n_seq << "\n";
    std::cout << "Number of links: " << info.graph_size.n_links << "\n";
    std::cout << "Total length: " << info.graph_size.total_length << "\n";
    std::cout << "Total length no overlaps: " << info.graph_size.total_length_no_overlaps << "\n\n";

    std::cout << "Sequences length info\n";
    for (const auto& perc : info.length_info.nx) {
        std::cout << "N" << perc.first << " : " << perc.second << "\n";
    }
    std::cout << "Average sequence length: " << info.length_info.average_length << "\n";
    std::cout << "Longest sequence: " << info.length_info.longest_sequence << "\n";
    std::cout << "Shortest sequence: " << info.length_info.shortest_sequence << "\n";
    for (const auto& perc : info.length_info.percentile_length) {
        std::cout << perc.first << " percentile of sequences length : " << perc.second << "\n";
    }
    std::cout << "\n";

    std::cout << "Graph connectivity info\n";
    std::cout << "Number of dead ends: " << info.connectivity_info.n_dead_ends << "\n";
    std::cout << "Percentage of dead ends: " << info.connectivity_info.dead_ends_percentage << "%\n";
    std::cout << "Number of connected components: " << info.connectivity_info.n_connected_components << "\n";
    std::cout << "Largest component: " << info.connectivity_info.largest_component_length << "\n";
    std::cout << "Max degree: " << info.connectivity_info.max_degree << "\n";
    std::cout << "Average degree: " << info.connectivity_info.average_degree << "\n\n";

    std::cout << "Coverage info\n";
    for (const auto& perc: info.coverage_info.percentile_coverage) {
        std::cout << perc.first << " percentile sequences coverage : " << perc.second << "\n";
    }
    std::cout << "Estimated genome length " << info.coverage_info.estimated_length << "\n";
}


void CalculateStat(const debruijn_graph::Graph& g, const std::vector<double>& n_percentiles,
               const std::vector<double>& median_length_percentiles,
               const std::vector<double>& sequences_coverage_percentiles,
               const std::filesystem::path& yaml_output_path) {
    StatInfo info;
    info.graph_size.k = g.k();

    size_t n_seq = 0;
    GetSequencesNumber(g, n_seq, info);

    GetLinksNumber(g, info);

    GetSequencesLengthInfo(g, median_length_percentiles, n_seq, info);

    GetDeadEndsInfo(g, info);

    GetConnectedComponentsInfo(g, info);

    for (size_t perc_ : n_percentiles) {
        info.length_info.nx[perc_] = Nx(g, perc_);
    }

    double median_coverage = 0;
    GetPercentileCoverage(g, sequences_coverage_percentiles, median_coverage, info);

    GetEstimatedGenomeLength(g, median_coverage, info);

    GetSequencesDegreeInfo(g, n_seq, info);

    if (!yaml_output_path.empty()) {
        if (yaml_output_path.extension() == ".yaml") {
            std::filesystem::path yaml_output = yaml_output_path;
            WriteToYaml(info, yaml_output);
        }
        else {
            WriteToYaml(info, std::cout);
        }
    }
    else {
        PrintStat(info);
    }
}