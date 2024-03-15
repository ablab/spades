//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <string>
#include <cmath>
#include <tuple>

#include "assembly_graph/core/graph.hpp"
#include "superpath_index.hpp"
#include "utils.hpp"

using debruijn_graph::EdgeId;
using debruijn_graph::ConjugateDeBruijnGraph;

class HMMPathInfo {
private:
    double rounded_score() const {
        return round(score * 1000.0) / 1000.0;
    }

    auto key() const {
        return std::make_tuple(rounded_score(), nuc_seq, label, path);
    }

    bool supported_by_original_path(const superpath_index::SuperpathIndex<EdgeId> &paths) const {
        if (label.size() == 0) {
            return false;
        }

        const std::string edge_prefix = "edge_";
        if (label.substr(0, edge_prefix.size()) == edge_prefix) {
            return false;
        }
        size_t path_id = std::stoll(label);
        size_t pos = superpath_index::find_subpath(path, paths[path_id]);
        return pos != size_t(-1);
    }

    friend void unique_hmm_path_info(std::vector<HMMPathInfo> &infos, const superpath_index::SuperpathIndex<EdgeId> &paths);

public:
    std::string hmmname;
    double score;
    std::string seq;
    std::string nuc_seq;
    std::vector<EdgeId> path;
    std::string alignment;
    std::string label;
    size_t pos;

    bool operator<(const HMMPathInfo &that) const {
        return key() < that.key();
    }

    size_t trim_first_edges(const ConjugateDeBruijnGraph &graph) {
        size_t edges_to_trim = 0;
        for (;edges_to_trim < path.size() - 1; ++edges_to_trim) {  // We cannot remove the last edge
            size_t edge_length_before_overlap = graph.length(path[edges_to_trim]);
            if (pos < edge_length_before_overlap) {
                break;
            }
            pos -= edge_length_before_overlap;
        }

        path.erase(path.begin(), path.begin() + edges_to_trim);
        return edges_to_trim;
    }

    HMMPathInfo(std::string name, double sc, std::string s, std::string nuc_s, std::vector<EdgeId> p, std::string alignment,
                std::string label, size_t pos)
            : hmmname(std::move(name)), score(sc),
              seq(std::move(s)), nuc_seq{std::move(nuc_s)},
              path(std::move(p)), alignment(std::move(alignment)), label{std::move(label)}, pos{pos} {}
};

inline void unique_hmm_path_info(std::vector<HMMPathInfo> &infos, const superpath_index::SuperpathIndex<EdgeId> &paths) {
    auto key = [&paths](const auto &info) {
        int supported = info.supported_by_original_path(paths) ? 0 : 1;
        return std::make_tuple(-info.rounded_score(), info.nuc_seq, info.path, supported, info.label);
    };
    sort_by(infos.begin(), infos.end(), key);

    auto eq = [](const auto &i1, const auto &i2) {
        return i1.rounded_score() == i2.rounded_score() && i1.nuc_seq == i2.nuc_seq && i1.path == i2.path;
    };

    auto it = std::unique(infos.begin(), infos.end(), eq);
    infos.erase(it, infos.end());
}
