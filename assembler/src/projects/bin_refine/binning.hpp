//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "pipeline/graph_pack.hpp"
#include <unordered_set>
#include <unordered_map>
#include <set>

namespace debruijn_graph {
class GraphPack;
}

namespace bin_stats {

class BinStats {
    using bin_id_t = std::string;
    using scaffold_id_t = std::string;
    using edge_binning_t = std::set<bin_id_t>;


    void ScaffoldsToEdges(const std::string& scaffolds_file,
                          const debruijn_graph::GraphPack &gp);

    static const std::string SCAFFOLD_NAME_PREFIX;
    static const std::string UNBINNED_ID;

 public:
    explicit BinStats(const debruijn_graph::Graph& g)
            : graph_(g) {}

    /// binning file in .tsv format (NODE_{scaffold_id}_* -> bin_id); scaffolds_file in .fasta format
    /// FIXME: Scaffolds should not be necessary
    void LoadBinning(const std::string& binning_file, const std::string& scaffolds_file,
                     const debruijn_graph::GraphPack &gp);

    const debruijn_graph::Graph& graph() const { return graph_;  }

    const auto& edges_binning() const { return edges_binning_; }
    auto& edges_binning() { return edges_binning_; }

    const auto& unbinned_edges() const { return unbinned_edges_; }
    auto& unbinned_edges() { return unbinned_edges_; }

    const auto& bins() const { return bins_; }

    friend std::ostream &operator<<(std::ostream &os, const BinStats &stats);
private:
    const debruijn_graph::Graph& graph_;
    std::unordered_map<scaffold_id_t, bin_id_t> scaffolds_binning_{};
    std::unordered_map<debruijn_graph::EdgeId, edge_binning_t> edges_binning_{};
    std::unordered_set<debruijn_graph::EdgeId> unbinned_edges_{};
    std::unordered_set<bin_id_t> bins_{};
};
}
