#pragma once

#include "aligner_output_reader.hpp"
#include "common.hpp"

namespace helpers {

inline double GetCov(std::string const & edge_name) {
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

inline debruijn_graph::EdgeId GetEdgeId(std::string const & edge_name, debruijn_graph::Graph const & graph) {
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

template<class Columns, Columns ... columns>
double GetIDY(Record<Columns, columns ...> const & record) {
    auto matched = record.template Get<Columns::match>();
    auto contig_start = record.template Get<Columns::Q_start>();
    auto contig_end = record.template Get<Columns::Q_end>();
    auto contig_delta = contig_end - contig_start;
    return static_cast<double>(matched) / static_cast<double>(contig_delta);
}

template<class Columns, Columns ... columns>
FilterType<Columns, columns ...> GetFilter(std::function<bool(debruijn_graph::EdgeId)> unique_edge_checker, debruijn_graph::Graph const & graph) {
    return [unique_edge_checker_ = std::move(unique_edge_checker), &graph](auto const & element) {
        constexpr auto EDGE_LENGTH_ERROR_COEFF = 0.05;
        constexpr auto LENGTHS_ERROR_COEFF = 0.05;
        constexpr auto MIN_TAIL_OVERLAP_LEN = 1000;
        const auto TAIL_GAP = static_cast<int>(graph.k());
        if (element.template Get<Columns::strand>() != '+')
            return false;
        auto contig_start = element.template Get<Columns::Q_start>();
        auto contig_end = element.template Get<Columns::Q_end>();
        auto contig_len = element.template Get<Columns::Q_size>();

        auto edge_start = element.template Get<Columns::T_start>();
        auto edge_end = element.template Get<Columns::T_end>();
        auto edge_len = element.template Get<Columns::T_size>();
        auto contig_delta = contig_end - contig_start;
        auto edge_delta = edge_end - edge_start;

        auto edge_len_difference = std::abs(edge_delta - edge_len);
        auto contig_len_difference = std::abs(contig_delta - contig_len);
        auto lens_difference = std::abs(contig_delta - edge_delta);

        auto const & edge_title = element.template Get<Columns::T_name>();
        
        bool edge_fully_inside_contig = (double)edge_len_difference <= (double)edge_len * EDGE_LENGTH_ERROR_COEFF;
        bool contig_fully_inside_edge = (double)contig_len_difference <= (double)contig_len * EDGE_LENGTH_ERROR_COEFF;

        bool contig_is_started_with_edge =
                MIN_TAIL_OVERLAP_LEN <= edge_delta &&
                edge_len - edge_end < TAIL_GAP &&
                contig_start < TAIL_GAP;

        bool contig_is_ended_with_edge =
                MIN_TAIL_OVERLAP_LEN <= edge_delta &&
                edge_start < TAIL_GAP &&
                contig_len - contig_end < TAIL_GAP;

        // if ((edge_title == "EDGE_557860_length_2034_cov_136.129864:EDGE_45218246_length_161_cov_375.679245;" || 
        //      edge_title == "EDGE_390252_length_1914_cov_124.899946':EDGE_503264_length_70_cov_461.800000;") &&
        //      element.template Get<Columns::Q_name>() == "contig_4")
        // {
        //     std::cout << element.template Get<Columns::Q_name>() << " -> " << element.template Get<Columns::T_name>() << ":\n";
        //     std::cout << "  " << unique_edge_checker_(GetEdgeId(edge_title, graph)) << '\n';
        //     std::cout << "  " << (GetIDY(element) > 0.90) << '\n';
        //     std::cout << "  " << ((double) lens_difference < (double) edge_delta * LENGTHS_ERROR_COEFF) << '\n';
        //     std::cout << "  " << (edge_fully_inside_contig || contig_fully_inside_edge || contig_is_started_with_edge || contig_is_ended_with_edge) << '\n';
        //     std::cout << "  " << (GetCov(edge_title) > 2) << '\n';
        //     return true;
        // }


        return unique_edge_checker_(GetEdgeId(edge_title, graph)) &&
               GetIDY(element) > 0.90 &&
               (double) lens_difference < (double) edge_delta * LENGTHS_ERROR_COEFF &&
               (edge_fully_inside_contig || contig_fully_inside_edge || contig_is_started_with_edge || contig_is_ended_with_edge) &&
               GetCov(edge_title) > 2;
    };
}

template<class Columns, Columns ... columns>
bool HasBadIntersection(Record<Columns, columns ...> const & lhs, Record<Columns, columns ...> const & rhs, size_t k) {
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

template<class Columns, Columns ... columns>
using DropAlg = std::function<size_t(std::vector<size_t> &, Records<Columns, columns ...> const &, size_t)>;

template<class Columns, Columns ... columns>
using DropMarker = std::function<std::pair<bool, bool>(Record<Columns, columns ...> const &, Record<Columns, columns ...> const &)>;

template<class Columns, Columns ... columns>
size_t DropperCore(std::vector<size_t> & contig, Records<Columns, columns ...> const & records, size_t k, DropMarker<Columns, columns ...> const & marker) {
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

template<class Columns, Columns ... columns>
DropAlg<Columns, columns ...> GetNonDropper() {
    return [](std::vector<size_t> &, Records<Columns, columns ...> const &, size_t) {
        return 0;
    };
}

template<class Columns, Columns ... columns>
DropAlg<Columns, columns ...> GetFullDropper() {
    return [](std::vector<size_t> & contig, Records<Columns, columns ...> const & records, size_t k) {
        DropMarker<Columns, columns ...> marker = [](auto const &, auto const &) { return std::make_pair(true, true); };
        return DropperCore(contig, records, k, marker);
    };
}

template<class Columns, Columns ... columns>
DropAlg<Columns, columns ...> GetLengthDropper() {
    return [](std::vector<size_t> & contig, Records<Columns, columns ...> const & records, size_t k) {
        DropMarker<Columns, columns ...> marker = [](auto const &lhs, auto const &rhs) {
            auto lhs_len = lhs.template Get<Columns::T_size>();
            auto rhs_len = rhs.template Get<Columns::T_size>();
            if (lhs_len >= 2*rhs_len)
                return std::make_pair(false, true);
            if (rhs_len >= 2*lhs_len)
                return std::make_pair(true, false);
            return std::make_pair(true, true);
        };
        return DropperCore(contig, records, k, marker);
    };
}

template<class Columns, Columns ... columns>
DropAlg<Columns, columns ...> GetTransitiveDropperByIDY() {
    return [](std::vector<size_t> & contig, Records<Columns, columns ...> const & records, size_t k) {
        DropMarker<Columns, columns ...> marker = [](auto const & lhs, auto const & rhs) {
            auto idy = GetIDY(lhs) < GetIDY(rhs);
            return std::make_pair(idy, !idy);
        };
        return DropperCore(contig, records, k, marker);
    };
}

template<class Columns, Columns ... columns>
MapFromContigNameToContigFragments GetContigFragments(Records<Columns, columns ...> const & records, DropAlg<Columns, columns ...> const & drop_alg, size_t k) {
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

template<class Columns, Columns ... columns>
PathWithEdgePostionsContainer MakePaths(Records<Columns, columns ...> const & records, MapFromContigNameToContigFragments const & contig_fragments, debruijn_graph::Graph const & graph) {
    PathWithEdgePostionsContainer paths;
    std::vector<std::string> paths_names;
    for (auto const & contig : contig_fragments) {
        PathWithEdgePostions path;
        for (size_t index : contig.second) {
            auto const & record = records[index];
            auto start_shift = record.template Get<Columns::T_start>() - 0;
            auto end_shift = record.template Get<Columns::T_size>() - record.template Get<Columns::T_end>();
            path.start_positions.push_back(record.template Get<Columns::Q_start>() - start_shift);
            path.end_positions.push_back(record.template Get<Columns::Q_end>() + end_shift);
            path.edges.push_back(GetEdgeId(record.template Get<Columns::T_name>(), graph));
        }
        if (!path.edges.empty()) {
            path.path_name = contig.first;
            paths.push_back(std::move(path));
        }
    }
    return paths;
}

} // namespace helpers
