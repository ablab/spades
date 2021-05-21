//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "binning.hpp"

#include "io/utils/id_mapper.hpp"

#include "io/reads/io_helper.hpp"
#include "utils/filesystem/path_helper.hpp"

using namespace debruijn_graph;
using namespace bin_stats;

const std::string Binning::UNBINNED_ID = "0";

EdgeLabels::EdgeLabels(const EdgeId e, const Binning& bin_stats)
        : e(e) {
    auto bins = bin_stats.edges_binning().find(e);
    is_binned = bins != bin_stats.edges_binning().end();
    labels_probabilities.resize(bin_stats.bins().size());

    is_repetitive = false;
    if (is_binned) {
        size_t sz = bins->second.size();
        //size_t sz = bin_stats.multiplicities().at(e);
        for (bin_stats::Binning::BinId bin : bins->second)
            labels_probabilities.set(bin, 1.0 / static_cast<double>(sz));
        is_repetitive = bin_stats.multiplicities().at(e) > 1;
    }
}

void Binning::ScaffoldsToEdges() {
  edges_binning_.clear();
  unbinned_edges_.clear();

  for (const auto &scaffold_entry : scaffolds_binning_) {
      ScaffoldId scaffold = scaffold_entry.first;
      BinId bin_id = scaffold_entry.second;

      for (EdgeId e : scaffolds_paths_.at(scaffold)) {
          if (bin_id != UNBINNED) {
              edges_binning_[e].insert(bin_id);
              edges_binning_[graph_.conjugate(e)].insert(bin_id);
          }
          edges_multiplicity_[e] += 1;
          edges_multiplicity_[graph_.conjugate(e)] += 1;
      }
  }

  // find unbinned edges
  for (EdgeId edge : graph_.edges()) {
    if (edges_binning_.count(edge))
      continue;

    unbinned_edges_.insert(edge);
    unbinned_edges_.insert(graph_.conjugate(edge));
  }
}

void Binning::LoadBinning(const std::string &binning_file) {
  fs::CheckFileExistenceFATAL(binning_file);

  scaffolds_binning_.clear();
  bins_.clear();
  std::ifstream binning_reader(binning_file);

  BinId max_bin_id = 0;
  for (std::string line; std::getline(binning_reader, line, '\n');) {
    std::string scaffold_name;
    BinLabel bin_label;

    std::istringstream line_stream(line);
    line_stream >> scaffold_name;
    line_stream >> bin_label;
    BinId cbin_id;
    if (bin_label == UNBINNED_ID) // unbinned scaffold
        cbin_id = UNBINNED;
    else {
        auto entry = bins_.find(bin_label);
        if (entry == bins_.end()) { // new bin label
            cbin_id = max_bin_id++;
            bin_labels_.emplace(cbin_id, bin_label);
            bins_.emplace(bin_label, cbin_id);
        } else {
            cbin_id = entry->second;
        }
    }

    auto scaffold_entry = scaffolds_.find(scaffold_name);
    if (scaffold_entry == scaffolds_.end()) {
        INFO("Unknown scaffold: " << scaffold_name);
        continue;
    }
    
    scaffolds_binning_[scaffold_entry->second] = cbin_id;
  }

  ScaffoldsToEdges();
}

void Binning::WriteToBinningFile(const std::string& binning_file,
                                  const SoftBinsAssignment &soft_edge_labels, const BinningAssignmentStrategy& assignment_strategy,
                                  const io::IdMapper<std::string> &edge_mapper) {
    std::ofstream out_tsv(binning_file + ".tsv");
    std::ofstream out_lens(binning_file + ".bin_weights");
    std::ofstream out_edges(binning_file + ".edge_weights");

    auto weight_sorter = [] (const auto &lhs, const auto &rhs) {
        if (math::eq(rhs.second, lhs.second))
            return lhs.first < rhs.first;

        return rhs.second < lhs.second;
    };

    std::vector<std::pair<ScaffoldId, size_t>> scaffolds_by_length;
    for (const auto &path_entry : scaffolds_paths_) {
        size_t length = graph_.k();
        for (EdgeId e : path_entry.second)
            length += graph_.length(e);

        scaffolds_by_length.emplace_back(path_entry.first, length);
    }

    std::sort(scaffolds_by_length.begin(), scaffolds_by_length.end(),
              [](const auto &lhs, const auto &rhs) {
                  if (lhs.second == rhs.second)
                      return lhs.first < lhs.first;

                  return lhs.second > rhs.second;
              });

    for (const auto &entry : scaffolds_by_length) {
        const std::string &scaffold_name = scaffolds_labels_.at(entry.first);
        const auto &bins_weights = scaffolds_bin_weights_.at(entry.first);
        std::vector<BinId> new_bin_id =
                assignment_strategy.ChooseMajorBins(bins_weights,
                                                    soft_edge_labels, *this);
        out_tsv << scaffold_name;
        if (new_bin_id.empty())
            out_tsv << '\t' << UNBINNED_ID;
        else
            for (BinId bin : new_bin_id)
                out_tsv << '\t' << bin_labels_.at(bin);
        out_tsv << '\n';
        out_lens << scaffold_name << '\t';
        out_lens << "nz: " << bins_weights.nonZeros();
        std::vector<std::pair<BinId, double>> weights;
        for (const auto &entry : bins_weights)
            weights.emplace_back(entry.index(), entry.value());
        std::sort(weights.begin(), weights.end(), weight_sorter);
        for (const auto &entry : weights)
            out_lens << '\t' << bin_labels_.at(entry.first) << ":" << entry.second;

        out_lens << '\n';
    }

    out_edges.precision(3);
    out_edges << "# edge_id\tinternal_edge_id\tlength\tbinned\twas_binned\trepetitive\tedge probs\n";
    for (EdgeId e : graph_.canonical_edges()) {
        const EdgeLabels& edge_labels = soft_edge_labels.at(e);
        out_edges << edge_mapper[graph_.int_id(e)] << '\t' << e << '\t' << graph_.length(e) + graph_.k()
                  << !unbinned_edges_.count(e) << '\t' << edge_labels.is_binned << '\t' << edge_labels.is_repetitive << '\t'
                  << "nz: " << edge_labels.labels_probabilities.nonZeros();
        std::vector<std::pair<BinId, double>> weights;
        for (const auto &entry : edge_labels.labels_probabilities)
            weights.emplace_back(entry.index(), entry.value());
        std::sort(weights.begin(), weights.end(), weight_sorter);
        for (const auto &entry : weights)
            out_edges << '\t' << bin_labels_.at(entry.first) << ":" << entry.second;
        out_edges << '\n';
    }
}

void Binning::AssignBins(const SoftBinsAssignment& soft_edge_labels,
                         const BinningAssignmentStrategy& assignment_strategy) {
    assignment_strategy.AssignEdgeBins(soft_edge_labels, *this);
    for (const auto &path_entry : scaffolds_paths_) {
        std::vector<EdgeId> scaffold_path(path_entry.second.begin(), path_entry.second.end());
        scaffolds_bin_weights_[path_entry.first] =
                assignment_strategy.AssignScaffoldBins(scaffold_path,
                                                       soft_edge_labels, *this);
    }
}

namespace bin_stats {
std::ostream &operator<<(std::ostream &os, const Binning &stats) {
    const auto &graph = stats.graph();

    os << "Total bins: " << stats.bins().size() << std::endl;
    size_t sum_length = 0;
    size_t unbinned_length = 0;
    for (const EdgeId e : graph.edges()) {
      const size_t length = graph.length(e);
      sum_length += length;
      if (stats.unbinned_edges().count(e)) {
        unbinned_length += length;
      }
    }

    os << "Total edges: " << graph.e_size()
       << ", binned: " << stats.edges_binning().size()
       << ", unbinned: " << stats.unbinned_edges().size()
       << " (" << (static_cast<double>(stats.unbinned_edges().size()) / static_cast<double>(graph.e_size()) * 100) << "%)" << std::endl
       << "Sum edge length: " << sum_length << std::endl
       << "Unbinned edges sum length: " << unbinned_length
       << " (" << (static_cast<double>(unbinned_length) / static_cast<double>(sum_length) * 100) << "%)";

    return os;
}

std::ostream &operator<<(std::ostream &os, const EdgeLabels &labels) {
    os << labels.is_binned << '\t';
    os << "nz: " << labels.labels_probabilities.nonZeros();
    for (const auto &entry : labels.labels_probabilities)
        os << '\t' << entry.index() << ":" << entry.value();

    return os;
}

} // namespace bin_stats
