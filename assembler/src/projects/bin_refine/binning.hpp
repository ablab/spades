//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"

#include "blaze/math/CompressedVector.h"

#include "binning_assignment_strategy.hpp"

#include <unordered_set>
#include <unordered_map>
#include <set>

namespace io {
template<typename T>
class IdMapper;
}

namespace bin_stats {

class BinningAssignmentStrategy;

class BinStats;

using LabelProbabilities = blaze::CompressedVector<double>;

struct EdgeLabels {
    // TODO: Could pack e and is_binned into single 64 bits
    debruijn_graph::EdgeId e;
    bool is_binned;
    LabelProbabilities labels_probabilities;

    EdgeLabels(debruijn_graph::EdgeId e, const BinStats& bin_stats);
    EdgeLabels(const EdgeLabels& edge_labels) = default;
    EdgeLabels& operator=(const EdgeLabels& edge_labels) = default;

    friend std::ostream &operator<<(std::ostream &os, const EdgeLabels &labels);
};

using SoftBinsAssignment = std::unordered_map<debruijn_graph::EdgeId, EdgeLabels>;

class BinStats {
    static const std::string UNBINNED_ID;
 public:
    using BinLabel = std::string;
    using BinId = uint64_t;
    using ScaffoldName = std::string;
    using EdgeBinning = std::unordered_set<BinId>;
    using ScaffoldsPaths = std::unordered_map<ScaffoldName, std::unordered_set<debruijn_graph::EdgeId>>;
    using BinLabels = std::unordered_map<BinId, BinLabel>;
    
    static constexpr BinId UNBINNED = BinId(-1);

    explicit BinStats(const debruijn_graph::Graph& g)
            : graph_(g) {}

    /// binning file in .tsv format (NODE_{scaffold_id}_* -> bin_id); scaffolds_file in .fasta format
    void LoadBinning(const std::string& binning_file, const ScaffoldsPaths &scaffolds_paths);
    void WriteToBinningFile(const std::string& binning_file, const ScaffoldsPaths &scaffolds_paths,
                            const SoftBinsAssignment &edge_soft_labels, const BinningAssignmentStrategy& assignment_strategy,
                            const io::IdMapper<std::string> &edge_mapper);
    void AssignEdgeBins(const SoftBinsAssignment& soft_bins_assignment, const BinningAssignmentStrategy& assignment_strategy);

    const debruijn_graph::Graph& graph() const { return graph_;  }

    const auto& edges_binning() const { return edges_binning_; }
    auto& edges_binning() { return edges_binning_; }

    const auto& unbinned_edges() const { return unbinned_edges_; }
    auto& unbinned_edges() { return unbinned_edges_; }

    const auto& bins() const { return bins_; }
    const auto& bin_labels() const { return bin_labels_; }

    friend std::ostream &operator<<(std::ostream &os, const BinStats &stats);
 private:
    void ScaffoldsToEdges(const ScaffoldsPaths &scaffolds_paths);

    const debruijn_graph::Graph& graph_;

    std::unordered_map<ScaffoldName, BinId> scaffolds_binning_{};
    BinLabels bin_labels_{};
    std::unordered_map<BinLabel, BinId> bins_{};
    std::unordered_map<debruijn_graph::EdgeId, EdgeBinning> edges_binning_{};
    std::unordered_set<debruijn_graph::EdgeId> unbinned_edges_{};
};
}
