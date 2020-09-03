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
    pos = edge_name.find(':');
    if (pos == 0)
        throw "Bad edge format at: " + edge_name;
    bool should_be_reversed = ((pos == std::string::npos ? edge_name.back() : edge_name[pos - 1]) == '\'');
    if (should_be_reversed)
        return graph.conjugate(edge_id);
    return edge_id;
}

FilterType<Columns::match, Columns::strand, Columns::block_count,
           Columns::Q_name, Columns::Q_size, Columns::Q_start, Columns::Q_end,
           Columns::T_name, Columns::T_size, Columns::T_start, Columns::T_end
           > GetFilter()
{
    return [](auto const & element) {
        if (element.template Get<Columns::strand>() != '+' || element.template Get<Columns::block_count>() != 1)
            return false;
        auto matched = element.template Get<Columns::match>();
        auto contig_start = element.template Get<Columns::Q_start>();
        auto contig_end = element.template Get<Columns::Q_end>();
        auto contig_len = element.template Get<Columns::Q_size>();

        auto edge_start = element.template Get<Columns::T_start>();
        auto edge_end = element.template Get<Columns::T_end>();
        auto edge_len = element.template Get<Columns::T_size>();
        auto contig_delta = contig_end - contig_start;
        auto edge_delta = edge_end - edge_start;

        auto identity = static_cast<double>(matched) / static_cast<double>(contig_delta + 1);
        auto cov = static_cast<double>(edge_delta) / static_cast<double>(edge_len);

        auto contig_first_or_last = contig_start < 10 || (contig_len - contig_end) < 20;
        auto edge_first_or_last = edge_start < 20 || (edge_len - edge_end) < 20;

        return ((contig_first_or_last && edge_first_or_last) || (cov > 0.98 && identity > 0.99 && contig_delta > 999)) &&
                (double) contig_delta > (double) edge_delta * 0.99 &&
                GetCov(element.template Get<Columns::T_name>()) > 2;
    };
}

template<Columns ... columns>
std::pair<PathWithEdgePostionsContainer, std::vector<std::string>> MakePaths(Records<columns ...> const & records, debruijn_graph::Graph const & graph) {
    std::unordered_map<std::string, std::vector<size_t>> ContigFragments;
    for (size_t i = 0; i < records.size(); ++i)
        ContigFragments[records[i].template Get<Columns::Q_name>()].push_back(i);


    for (auto & contig : ContigFragments) {
        auto by_start_pos = [&records](size_t lhs, size_t rhs) {
            return records[lhs].template Get<Columns::Q_start>() < records[rhs].template Get<Columns::Q_start>();
        };
        std::sort(contig.second.begin(), contig.second.end(), by_start_pos);
    }

    PathWithEdgePostionsContainer paths;
    std::vector<std::string> paths_names;

    for (auto const & contig : ContigFragments) {
        PathWithEdgePostions path;
        for (size_t index : contig.second) {
            path.positions.push_back(records[index].template Get<Columns::Q_start>());
            path.edge_set.push_back({GetEdgeId(records[index].template Get<Columns::T_name>(), graph)});
        }
        paths.push_back(std::move(path));
        paths_names.push_back(contig.first);
    }
    return {std::move(paths), std::move(paths_names)};
}
