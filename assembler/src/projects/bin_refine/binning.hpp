//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "binning_assignment_strategy.hpp"
#include "id_map.hpp"

#include "assembly_graph/core/graph.hpp"

#include <blaze/math/CompressedVector.h>
#include <unordered_set>
#include <unordered_map>
#include <set>

namespace io {
template<typename T>
class IdMapper;
}

namespace bin_stats {

class BinningAssignmentStrategy;

class Binning;

using LabelProbabilities = blaze::CompressedVector<double>;

enum OutputOptions : uint64_t {
    CAMI        = 1 << 0,
    EmitZeroBin = 1 << 1,
    TallMulti   = 1 << 2,
};

struct EdgeLabels {
    using BinId = uint64_t;

    static constexpr BinId UNBINNED = BinId(-1);

    // TODO: Could pack e and is_binned into single 64 bits
    debruijn_graph::EdgeId e;
    bool is_binned : 1;
    bool is_repetitive : 1;
    LabelProbabilities labels_probabilities;

    EdgeLabels()
            : e(0), is_binned(false), is_repetitive(false) {}
    EdgeLabels(const debruijn_graph::EdgeId &e,
               bool is_binned,
               bool is_repetitive,
               const LabelProbabilities &labels_probabilities) : e(e),
                                                                 is_binned(is_binned),
                                                                 is_repetitive(is_repetitive),
                                                                 labels_probabilities(labels_probabilities) {}

    EdgeLabels(debruijn_graph::EdgeId e, const Binning& binning);
    EdgeLabels(const EdgeLabels& edge_labels) = default;
    EdgeLabels& operator=(const EdgeLabels& edge_labels) = default;

    friend std::ostream &operator<<(std::ostream &os, const EdgeLabels &labels);
};

struct BinStats {
    double mean_cov;
    double sd_cov;
    double m2_cov;
    double len;
};

using SoftBinsAssignment = adt::id_map<EdgeLabels, debruijn_graph::EdgeId>;

class Binning {
    static const std::string UNBINNED_ID;
 public:
    using BinLabel = std::string;
    using BinId = EdgeLabels::BinId;
    using ScaffoldName = std::string;
    using ScaffoldId = uint64_t;
    using ScaffoldPath = std::unordered_set<debruijn_graph::EdgeId>;
    using Scaffold = std::pair<ScaffoldName, ScaffoldPath>;
    using EdgeBinning = std::unordered_set<BinId>;
    using ScaffoldsPaths = std::unordered_map<ScaffoldId, ScaffoldPath>;
    using BinLabels = std::unordered_map<BinId, BinLabel>;
    using ScaffoldLabels = std::unordered_map<ScaffoldId, ScaffoldName>;

    static constexpr BinId UNBINNED = EdgeLabels::UNBINNED;

    explicit Binning(const debruijn_graph::Graph& g)
            : graph_(g) {}

    void InitScaffolds(const std::vector<Scaffold> &scaffold_paths) {
        for (ScaffoldId id = 0; id < scaffold_paths.size(); ++id) {
            scaffolds_.emplace(scaffold_paths[id].first, id);
            scaffolds_labels_.emplace(id, scaffold_paths[id].first);
            scaffolds_paths_.emplace(id, scaffold_paths[id].second);
        }
    }

    /// binning file in .tsv format (NODE_{scaffold_id}_* -> bin_id); scaffolds_file in .fasta format
    void LoadBinning(const std::string& binning_file, bool cami);
    void AssignBins(const SoftBinsAssignment& soft_bins_assignment, const BinningAssignmentStrategy& assignment_strategy);

    blaze::DynamicMatrix<double> BinDistance(const SoftBinsAssignment& soft_bins_assignment,
                                             bool edges = false);

    void WriteToBinningFile(const std::string& binning_file, uint64_t output_options,
                            const SoftBinsAssignment &edge_soft_labels, const BinningAssignmentStrategy& assignment_strategy,
                            const io::IdMapper<std::string> &edge_mapper);

    const debruijn_graph::Graph& graph() const { return graph_;  }

    const auto& edges_binning() const { return edges_binning_; }
    auto& edges_binning() { return edges_binning_; }

    const auto& unbinned_edges() const { return unbinned_edges_; }
    auto& unbinned_edges() { return unbinned_edges_; }

    const auto& bins() const { return bins_; }
    const auto& bin_labels() const { return bin_labels_; }

    const auto& multiplicities() const { return edges_multiplicity_; }

    friend std::ostream &operator<<(std::ostream &os, const Binning &binning);
 private:
    void ScaffoldsToEdges();

    const debruijn_graph::Graph& graph_;
    std::string sample_id_ = "";

    // All about scaffolds
    std::unordered_map<ScaffoldId, BinId> scaffolds_binning_{};
    std::unordered_map<ScaffoldId, LabelProbabilities> scaffolds_bin_weights_{};
    ScaffoldsPaths scaffolds_paths_{};
    ScaffoldLabels scaffolds_labels_{};
    std::unordered_map<ScaffoldName, ScaffoldId> scaffolds_{};

    // All about bins
    std::unordered_map<BinLabel, BinId> bins_{};
    BinLabels bin_labels_{};
    std::unordered_map<BinId, BinStats> bin_stats_{};

    // All about edges
    // FIXME: id_map
    std::unordered_map<debruijn_graph::EdgeId, EdgeBinning> edges_binning_{};
    std::unordered_set<debruijn_graph::EdgeId> unbinned_edges_{};
    std::unordered_map<debruijn_graph::EdgeId, size_t> edges_multiplicity_{};
};

class LabelInitializer {
  public:
    using Graph = debruijn_graph::Graph;

    LabelInitializer(const Graph &g);

    SoftBinsAssignment InitLabels(const Binning &bin_stats) const;
  private:
    const Graph &g_;
};
}
