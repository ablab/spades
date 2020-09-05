#pragma once

#include "blat_output_reader.hpp"
#include "common.hpp"

double GetCov(std::string const & edge_name) {
    std::string cov_pattern = "_cov_";
    auto pos = edge_name.find(cov_pattern);
    if (pos == std::string::npos)
        throw "No cov info at: " + edge_name;

    pos += cov_pattern.size();
    size_t len = 0;
    while (pos + len < edge_name.size() && isdigit(edge_name[pos + len]))
        ++len;
    if (pos + len + 1 < edge_name.size() && edge_name[pos + len] == '.' && isdigit(edge_name[pos + len + 1])) {
        len += 2;
        while (pos + len < edge_name.size() && isdigit(edge_name[pos + len]))
            ++len;
    }
    if (len == 0)
        throw "Bad cov format at: " + edge_name;
    return std::stod(edge_name.substr(pos, len));
}

debruijn_graph::EdgeId GetEdgeId(std::string const & edge_name, debruijn_graph::Graph const & graph) {
    std::string cov_pattern = "EDGE_";
    auto pos = edge_name.find(cov_pattern);
    if (pos == std::string::npos)
        throw "No edge id at: " + edge_name;

    pos += cov_pattern.size();
    size_t len = 0;
    while (pos + len < edge_name.size() && isdigit(edge_name[pos + len]))
        ++len;
    if (len == 0)
        throw "Bad edge format at: " + edge_name;
    debruijn_graph::EdgeId edge_id = std::stoull(edge_name.substr(pos, len));
    pos = edge_name.find_first_of(":;");
    if (pos == 0 || pos == std::string::npos)
        throw "Bad edge format at: " + edge_name;
    if (edge_name[pos - 1] == '\'')
        return graph.conjugate(edge_id);
    return edge_id;
}

template<Columns ... columns>
double GetIDY(Record<columns ...> const & record) {
    auto matched = record.template Get<Columns::match>();
    auto contig_start = record.template Get<Columns::Q_start>();
    auto contig_end = record.template Get<Columns::Q_end>();
    auto contig_delta = contig_end - contig_start + 1;
    return static_cast<double>(matched) / static_cast<double>(contig_delta);
}

template<Columns ... columns>
FilterType<columns ...> GetFilter() {
    return [](auto const & element) {
        constexpr auto EDGE_LENGTH_ERROR_COEFF = 0.02;
        constexpr auto LENGTHS_ERROR_COEFF = 0.01;
        if (element.template Get<Columns::strand>() != '+' || element.template Get<Columns::block_count>() != 1)
            return false;
        auto contig_start = element.template Get<Columns::Q_start>();
        auto contig_end = element.template Get<Columns::Q_end>();

        auto edge_start = element.template Get<Columns::T_start>();
        auto edge_end = element.template Get<Columns::T_end>();
        auto edge_len = element.template Get<Columns::T_size>();
        auto contig_delta = contig_end - contig_start + 1;
        auto edge_delta = edge_end - edge_start + 1;

        auto edge_len_difference = std::abs(edge_delta - edge_len);
        auto lens_difference = std::abs(contig_delta - edge_delta);

        return GetIDY(element) > 0.95 &&
               contig_delta > 100 &&
               (double)edge_len_difference < (double)edge_len * EDGE_LENGTH_ERROR_COEFF &&
               (double) lens_difference < (double) edge_delta * LENGTHS_ERROR_COEFF &&
               GetCov(element.template Get<Columns::T_name>()) > 2;
    };
}

template<Columns ... columns>
bool HasBadIntersection(Record<columns ...> const & lhs, Record<columns ...> const & rhs, size_t k) {
    auto lhs_start = lhs.template Get<Columns::Q_start>();
    auto lhs_end = lhs.template Get<Columns::Q_end>() + 1;
    auto rhs_start = rhs.template Get<Columns::Q_start>();
    auto rhs_end = rhs.template Get<Columns::Q_end>() + 1;
    if (rhs_start < lhs_start) {
        std::swap(lhs_start, rhs_start);
        std::swap(lhs_end, rhs_end);
    }
    k += 5;
    // now lhs_start < rhs_start
    return rhs_start + static_cast<long long>(k) < lhs_end || rhs_end <= lhs_end;
}

using MapFromContigNameToContigFragments = std::unordered_map<std::string, std::vector<size_t>>;

template<Columns ... columns>
using DropAlg = std::function<size_t(std::vector<size_t> &, Records<columns ...> const &, size_t)>;

template<Columns ... columns>
using DropMarker = std::function<std::pair<bool, bool>(Record<columns ...> const &, Record<columns ...> const &)>;

template<Columns ... columns>
size_t DropperCore(std::vector<size_t> & contig, Records<columns ...> const & records, size_t k, DropMarker<columns ...> const & marker) {
    std::vector<bool> is_bad_edge(contig.size(), false);
    for (size_t i = 0; i + 1 < contig.size(); ++i) {
        for (size_t j = i + 1; j < contig.size(); ++j) {
            if (HasBadIntersection(records[contig[i]], records[contig[j]], k)) {
                auto marks = marker(records[contig[i]], records[contig[j]]);
                is_bad_edge[i] = is_bad_edge[i] | marks.first;
                is_bad_edge[j] = is_bad_edge[j] | marks.second;
            }
        }
    }
    std::vector<size_t> new_fragments;
    for (size_t i = 0; i < contig.size(); ++i) {
        if (!is_bad_edge[i])
            new_fragments.push_back(contig[i]);
    }
    auto dropped_pats_cnt = contig.size() - new_fragments.size();
    contig = std::move(new_fragments);
    return dropped_pats_cnt;
}

template<Columns ... columns>
DropAlg<columns ...> GetNonDropper() {
    return [](std::vector<size_t> &, Records<columns ...> const &, size_t) {
        return 0;
    };
}

template<Columns ... columns>
DropAlg<columns ...> GetFullDropper() {
    return [](std::vector<size_t> & contig, Records<columns ...> const & records, size_t k) {
        DropMarker<columns ...> marker = [](auto const &, auto const &) { return std::make_pair(true, true); };
        return DropperCore(contig, records, k, marker);
    };
}

template<Columns ... columns>
DropAlg<columns ...> GetTransitiveDropperByIDY() {
    return [](std::vector<size_t> & contig, Records<columns ...> const & records, size_t k) {
        DropMarker<columns ...> marker = [](auto const & lhs, auto const & rhs) {
            auto idy = GetIDY(lhs) < GetIDY(rhs);
            return std::make_pair(idy, !idy);
        };
        return DropperCore(contig, records, k, marker);
    };
}

template<Columns ... columns>
MapFromContigNameToContigFragments GetContigFragments(Records<columns ...> const & records, DropAlg<columns ...> const & drop_alg, size_t k) {
    MapFromContigNameToContigFragments contig_fragments;
    for (size_t i = 0; i < records.size(); ++i)
        contig_fragments[records[i].template Get<Columns::Q_name>()].push_back(i);

    for (auto & contig : contig_fragments) {
        auto by_start_pos = [&records](size_t lhs, size_t rhs) {
            return records[lhs].template Get<Columns::Q_start>() < records[rhs].template Get<Columns::Q_start>();
        };
        std::sort(contig.second.begin(), contig.second.end(), by_start_pos);
        auto dropped = drop_alg(contig.second, records, k);
        if (dropped) 
            WARN(contig.first << ": " << dropped << " edge" << (dropped != 1 ? "s" : "") << " would be dropped");
    }
    return contig_fragments;
}

template<Columns ... columns>
std::pair<PathWithEdgePostionsContainer, std::vector<std::string>> MakePaths(Records<columns ...> const & records, MapFromContigNameToContigFragments const & contig_fragments, debruijn_graph::Graph const & graph) {
    PathWithEdgePostionsContainer paths;
    std::vector<std::string> paths_names;
    for (auto const & contig : contig_fragments) {
        PathWithEdgePostions path;
        for (size_t index : contig.second) {
            path.positions.push_back((int)records[index].template Get<Columns::Q_start>());
            path.edge_set.push_back({GetEdgeId(records[index].template Get<Columns::T_name>(), graph)});
        }
        if (!path.positions.empty()) {
            paths.push_back(std::move(path));
            paths_names.push_back(contig.first);
        }
    }
    return {std::move(paths), std::move(paths_names)};
}
