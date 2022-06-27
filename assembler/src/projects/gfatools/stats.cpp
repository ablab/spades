//***************************************************************************
//* Copyright (c) 2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "stats.hpp"

#include "common/adt/id_map.hpp"
#include "common/assembly_graph/components/connected_component.hpp"
#include "common/toolchain/subgraph_utils.hpp"

#include "llvm/Support/YAMLTraits.h"

namespace gfa_tools {

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
        std::map<double, uint64_t> percentile_length;
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
        double median_coverage;
        uint64_t estimated_length;
    };
    GraphSize graph_size;
    SequencesLengthInfo length_info;
    GraphConnectivityInfo connectivity_info;
    CoverageInfo coverage_info;
};
}

namespace llvm { namespace yaml {

template<> struct MappingTraits<gfa_tools::StatInfo> {
    static void mapping(IO& io, gfa_tools::StatInfo& info) {
        io.mapRequired("Graph size", info.graph_size);
        io.mapRequired("Sequences length info", info.length_info);
        io.mapRequired("Graph connectivity info", info.connectivity_info);
        io.mapRequired("Coverage info", info.coverage_info);
    }
};

template<> struct MappingTraits<gfa_tools::StatInfo::GraphSize> {
    static void mapping(IO& io, gfa_tools::StatInfo::GraphSize& info) {
        io.mapRequired("Sequences overlap", info.k);
        io.mapRequired("Number of sequences", info.n_seq);
        io.mapRequired("Number of links", info.n_links);
        io.mapRequired("Total length", info.total_length);
        io.mapRequired("Total length no overlaps", info.total_length_no_overlaps);
    }
};

template<> struct MappingTraits<gfa_tools::StatInfo::SequencesLengthInfo> {
    static void mapping(IO& io, gfa_tools::StatInfo::SequencesLengthInfo& info) {
        io.mapRequired("N stat for percentiles", info.nx);
        io.mapRequired("Average sequence length", info.average_length);
        io.mapRequired("Longest sequence", info.longest_sequence);
        io.mapRequired("Shortest sequence", info.shortest_sequence);
        io.mapRequired("Length for percentiles", info.percentile_length);
    }
};

template<> struct MappingTraits<gfa_tools::StatInfo::GraphConnectivityInfo> {
    static void mapping(IO& io, gfa_tools::StatInfo::GraphConnectivityInfo& info) {
        io.mapRequired("Number of dead ends", info.n_dead_ends);
        io.mapRequired("Percentage of dead ends", info.dead_ends_percentage);
        io.mapRequired("Number of connected components", info.n_connected_components);
        io.mapRequired("Largest component", info.largest_component_length);
        io.mapRequired("Max degree", info.max_degree);
        io.mapRequired("Average degree", info.average_degree);
    }
};

template<> struct MappingTraits<gfa_tools::StatInfo::CoverageInfo> {
    static void mapping(IO& io, gfa_tools::StatInfo::CoverageInfo& info) {
        io.mapRequired("Coverage for percentile", info.percentile_coverage);
        io.mapRequired("Median coverage", info.median_coverage);
        io.mapRequired("Estimated genome length", info.estimated_length);
    }
};
} }

namespace gfa_tools {

template<class T>
void WriteToYaml(T &t, const std::filesystem::path& filename) {
    std::ofstream ofs(filename, std::ios::binary);
    WriteToYaml(t, ofs);
    ofs.close();
}

template<class T>
void WriteToYaml(T &t, std::ostream& os) {
    std::string buffer;
    llvm::raw_string_ostream raw_os(buffer);
    llvm::yaml::Output yout(raw_os);
    yout << t;
    raw_os.str();
    os << buffer;
}

void CountLinksNumber(const debruijn_graph::Graph& g, StatInfo& info) {
    size_t n_links = 0;
    for (debruijn_graph::VertexId v : g.canonical_vertices()) {
        n_links += (g.IncomingEdgeCount(v) + g.OutgoingEdgeCount(v));
    }
    info.graph_size.n_links = n_links;
}


void CountDeadEndsInfo(const debruijn_graph::Graph& g, StatInfo& info) {
    size_t n_dead_ends = 0;
    for (debruijn_graph::EdgeId e : g.canonical_edges()) {
        if (toolchain::IsDeadEnd(g, g.EdgeStart(e)) || toolchain::IsDeadEnd(g, g.EdgeEnd(e))) {
            n_dead_ends += 1;
        }
    }
    info.connectivity_info.n_dead_ends = n_dead_ends;
    info.connectivity_info.dead_ends_percentage = n_dead_ends * 100.0 / g.e_size();
}


std::map<double, size_t> GetPercentile(const std::vector<size_t>& v, const std::vector<double>& percentiles) {
    VERIFY(std::is_sorted(v.begin(), v.end()));
    std::map<double, uint64_t> result;
    if (v.size() > 0) {
        for (double perc: percentiles) {
            if (perc < 0 || perc > 100) {
                continue;
            }
            size_t idx = floor((v.size() - 1) * perc / 100);
            result[perc] = v[idx];
        }
    }
    return result;
}


void CountSequencesLengthInfo(const debruijn_graph::Graph& g,
                              const std::vector<double>& percentiles, StatInfo& info) {
    std::vector<size_t> sequences_lengths;
    uint64_t total_length = 0;
    uint64_t total_length_no_overlaps = 0;
    for (debruijn_graph::EdgeId e: g.edges()) {
        size_t sequence_length = g.EdgeNucls(e).size();
        sequences_lengths.push_back(sequence_length);
        total_length += sequence_length;
    }

    /*
     * To calculate a total number of base pairs in the graph without overlaps, it is needed to
     * count overlaps. Every vertex has a number of incoming and outgoing edges. Let the sum of them
     * be vertex degree. Vertex degree is a number of repetitions of a k-mer in this vertex. In a
     * total length without overlaps should be one repetition of this k-mer, so we delete the rest of them.
     * But in case when a vertex has a buckle, k-mer shouldn't be subtracted from one edge twice.
     */
    uint64_t doubles = 0;
    for (debruijn_graph::VertexId v : g.vertices()) {
        doubles += (g.IncomingEdgeCount(v) + g.OutgoingEdgeCount(v) - 1);
    }
    for (debruijn_graph::EdgeId e : g.edges()) {
        if (g.EdgeStart(e) == g.EdgeEnd(e)) {
            doubles -= 1;
        }
    }
    /*
     * In the total length all edges are taken into account, including those for which
     * there are reverse-complementary edges. Therefore, the resulting length must be divided by 2.
     */
    info.graph_size.total_length = total_length / 2;
    info.graph_size.total_length_no_overlaps = (total_length - info.graph_size.k * doubles) / 2;

    info.length_info.average_length = 1.0 * total_length / (2 * info.graph_size.n_seq);

    std::sort(sequences_lengths.begin(), sequences_lengths.end());
    info.length_info.shortest_sequence = sequences_lengths[0];
    info.length_info.percentile_length = GetPercentile(sequences_lengths, percentiles);
    info.length_info.longest_sequence = sequences_lengths.back();
}


void CountConnectedComponentsInfo(const debruijn_graph::Graph& g, StatInfo& info) {
    debruijn_graph::ConnectedComponentCounter connected_component_counter(g);
    connected_component_counter.CalculateComponents();

    size_t n_components = 0;
    for (const auto& edge : connected_component_counter.component_ids_) {
        n_components = std::max(n_components, edge.second + 1);
    }
    info.connectivity_info.n_connected_components = n_components;

    std::unordered_map<size_t, size_t> components_length;
    for (debruijn_graph::EdgeId e : g.canonical_edges()) {
        components_length[connected_component_counter.component_ids_[e]] += g.EdgeNucls(e).size();
    }

    size_t largest_component_size = 0;
    for (const auto& component : components_length) {
        largest_component_size = std::max(largest_component_size, component.second);
    }
    info.connectivity_info.largest_component_length = largest_component_size;
}


size_t GetNx(const debruijn_graph::Graph& g, double percent) {
    size_t sum_edge_length = 0;
    std::vector<size_t> lengths;
    for (debruijn_graph::EdgeId e : g.canonical_edges()) {
        lengths.push_back(g.EdgeNucls(e).size());
        sum_edge_length += g.EdgeNucls(e).size();
    }
    std::sort(lengths.rbegin(), lengths.rend());
    double len_perc = percent * sum_edge_length / 100.0;
    for (size_t i = 0; i < lengths.size(); i++) {
        if (double(lengths[i]) >= len_perc) {
            return lengths[i];
        } else {
            len_perc -= (double) lengths[i];
        }
    }
    return 0;
}


void CountPercentileCoverage(const debruijn_graph::Graph& g, const std::vector<double>& percentiles, StatInfo& info) {
    /*
     * This statistic is calculated in a such way because it is important to know such a value
     * that no less than a half of base pairs have the bigger coverage.
     */
    VERIFY(std::is_sorted(percentiles.begin(), percentiles.end()));
    std::vector<std::pair<double, size_t>> edges;
    uint64_t total_length = 0;
    for (debruijn_graph::EdgeId e: g.canonical_edges()) {
        edges.emplace_back(g.coverage(e), g.EdgeNucls(e).size());
        total_length += g.EdgeNucls(e).size();
    }
    std::sort(edges.begin(), edges.end());
    uint64_t cumulative_length = 0;
    size_t perc_ind = 0;
    bool median_is_counted = false;
    for (auto e : edges) {
        cumulative_length += e.second;
        while (perc_ind < percentiles.size() && cumulative_length >= total_length * percentiles[perc_ind] / 100.0) {
            info.coverage_info.percentile_coverage[percentiles[perc_ind++]] = e.first;
        }
        if (!median_is_counted && cumulative_length >= total_length * 0.5) {
            info.coverage_info.median_coverage = e.first;
            median_is_counted = true;
        }
    }
}


void CountEstimatedGenomeLength(const debruijn_graph::Graph& g, StatInfo& info) {
    uint64_t estimated_length = 0;
    adt::id_map<size_t, debruijn_graph::VertexId> vertexes_degrees(g.size());
    for (debruijn_graph::VertexId v : g.vertices()) {
        vertexes_degrees[v] = g.IncomingEdgeCount(v) + g.OutgoingEdgeCount(v);
    }
    for (debruijn_graph::EdgeId e: g.edges()) {
        size_t edge_length = g.EdgeNucls(e).size();
        if (edge_length > info.graph_size.k && vertexes_degrees[g.EdgeStart(e)] > 1) {
            edge_length -= info.graph_size.k;
            vertexes_degrees[g.EdgeStart(e)] -= 1;
        }
        if (g.EdgeStart(e) != g.EdgeEnd(e) && edge_length > info.graph_size.k && vertexes_degrees[g.EdgeEnd(e)] > 1) {
            edge_length -= info.graph_size.k;
            vertexes_degrees[g.EdgeEnd(e)] -= 1;
        }
        double coverage = g.coverage(e);
        size_t relative_coverage = floor(coverage / info.coverage_info.median_coverage);
        estimated_length += edge_length * relative_coverage;
    }
    info.coverage_info.estimated_length = estimated_length / 2;
}


void CountSequencesDegreeInfo(const debruijn_graph::Graph& g, StatInfo& info) {
    size_t max_degree = 0;
    double average_degree = 0;
    for (debruijn_graph::VertexId v : g.vertices()) {
        max_degree = std::max({max_degree, g.IncomingEdgeCount(v) + g.OutgoingEdgeCount(v)});
        average_degree += (g.IncomingEdgeCount(v) + g.OutgoingEdgeCount(v));
    }

    info.connectivity_info.max_degree = max_degree;
    info.connectivity_info.average_degree = average_degree / g.size();
}

void CalculateStat(const debruijn_graph::Graph& g, const std::vector<double>& n_percentiles,
                   const std::vector<double>& median_length_percentiles,
                   const std::vector<double>& sequences_coverage_percentiles,
                   const std::filesystem::path& yaml_output_path) {
    StatInfo info;
    info.graph_size.k = g.k();
    info.graph_size.n_seq = std::distance(g.e_begin<true>(), g.e_end<true>());

    CountLinksNumber(g, info);

    CountSequencesLengthInfo(g, median_length_percentiles, info);

    CountDeadEndsInfo(g, info);

    CountConnectedComponentsInfo(g, info);

    for (size_t perc : n_percentiles) {
        info.length_info.nx[perc] = GetNx(g, perc);
    }

    CountPercentileCoverage(g, sequences_coverage_percentiles, info);

    CountEstimatedGenomeLength(g, info);

    CountSequencesDegreeInfo(g, info);

    if (yaml_output_path.empty() || yaml_output_path == "-") {
        WriteToYaml(info, std::cout);
    } else {
        WriteToYaml(info, yaml_output_path);
        INFO("Information is written to " << absolute(yaml_output_path));
    }
}
}
