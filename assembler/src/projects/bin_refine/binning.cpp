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

const std::string BinStats::UNBINNED_ID = "0";

EdgeLabels::EdgeLabels(const EdgeId e, const BinStats& bin_stats)
        : e(e) {
    auto bins = bin_stats.edges_binning().find(e);
    is_binned = bins != bin_stats.edges_binning().end();
    labels_probabilities.resize(bin_stats.bins().size());

    if (is_binned) {
        size_t sz = bins->second.size();
        for (bin_stats::BinStats::BinId bin : bins->second)
            labels_probabilities.set(bin, 1.0 / static_cast<double>(sz));
    }
}

void BinStats::ScaffoldsToEdges(const ScaffoldsPaths &scaffolds_paths) {
  edges_binning_.clear();
  unbinned_edges_.clear();

  for (const auto &scaffold_entry : scaffolds_binning_) {
      const std::string& scaffold_name = scaffold_entry.first;

      BinId bin_id = scaffold_entry.second;
      auto path_entry = scaffolds_paths.find(scaffold_name);
      if (path_entry == scaffolds_paths.end()) {
          INFO("No path for scaffold " << scaffold_name);
          continue;
      }

      for (EdgeId e : path_entry->second) {
          edges_binning_[e].insert(bin_id);
          edges_binning_[graph_.conjugate(e)].insert(bin_id);
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

void BinStats::LoadBinning(const std::string& binning_file, const ScaffoldsPaths &scaffolds_paths) {
  fs::CheckFileExistenceFATAL(binning_file);

  scaffolds_binning_.clear();
  bins_.clear();
  std::ifstream binning_reader(binning_file);

  BinId max_bin_id = 0;
  for (std::string line; std::getline(binning_reader, line, '\n');) {
    std::string scaffold_name;
    std::istringstream line_stream(line);
    line_stream >> scaffold_name;
    BinLabel bin_label;
    line_stream >> bin_label;
    if (bin_label == UNBINNED_ID) // unbinned scaffold, skip
      continue;

    BinId cbin_id;
    auto entry = bins_.find(bin_label);
    if (entry == bins_.end()) { // new bin label
        cbin_id = max_bin_id++;
        bin_labels_.emplace(cbin_id, bin_label);
        bins_.emplace(bin_label, cbin_id);
    } else {
        cbin_id = entry->second;
    }

    scaffolds_binning_[scaffold_name] = cbin_id;
  }

  ScaffoldsToEdges(scaffolds_paths);
}

void BinStats::WriteToBinningFile(const std::string& binning_file, const ScaffoldsPaths &scaffolds_paths,
                                  const SoftBinsAssignment &soft_edge_labels, const BinningAssignmentStrategy& assignment_strategy,
                                  const io::IdMapper<std::string> &edge_mapper) {
    std::ofstream out_tsv(binning_file + ".tsv");
    std::ofstream out_lens(binning_file + ".bin_weights");
    std::ofstream out_edges(binning_file + ".edge_weights");

    for (const auto &path_entry : scaffolds_paths) {
      const std::string& scaffold_name = path_entry.first;

      std::vector<EdgeId> scaffold_path(path_entry.second.begin(), path_entry.second.end());
      auto bins_weights = assignment_strategy.AssignScaffoldBins(scaffold_path,
                                                                 soft_edge_labels, *this);
      std::vector<BinId> new_bin_id = assignment_strategy.ChooseMajorBins(bins_weights,
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
      for (const auto &entry : bins_weights)
          out_lens << '\t' << entry.index() << ":" << entry.value();
      out_lens << '\n';
    }

    out_edges.precision(3);
    out_edges << "# edge_id\tbinned\twas_binned\tedge probs\n";
    for (EdgeId e : graph_.canonical_edges()) {
        const EdgeLabels& edge_labels = soft_edge_labels.at(e);
        out_edges << edge_mapper[graph_.int_id(e)] << '\t' << !unbinned_edges_.count(e) << '\t' << edge_labels << '\n';
    }
}

void BinStats::AssignEdgeBins(const SoftBinsAssignment& soft_bins_assignment, const BinningAssignmentStrategy& assignment_strategy) {
    assignment_strategy.AssignEdgeBins(soft_bins_assignment, *this);
}

namespace bin_stats {
std::ostream &operator<<(std::ostream &os, const BinStats &stats) {
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

    os << "Total edges: " << graph.e_size() << std::endl
       << "Total edges length: " << sum_length << std::endl
       << "Total binned edges: " << stats.edges_binning().size() << std::endl
       << "Unbinned edges: " << stats.unbinned_edges().size() << std::endl
       << "Unbinned edges total length: " << unbinned_length << std::endl
       << "Unbinned edges number ratio: " << static_cast<double>(stats.unbinned_edges().size()) / static_cast<double>(graph.e_size()) << std::endl
       << "Unbinned edges length ratio: " << static_cast<double>(unbinned_length) / static_cast<double>(sum_length) << std::endl;

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