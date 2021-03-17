//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "binning.hpp"

#include "modules/alignment/sequence_mapper.hpp"

#include "io/reads/io_helper.hpp"
#include "utils/filesystem/path_helper.hpp"

using namespace debruijn_graph;
using namespace bin_stats;

const std::string BinStats::UNBINNED_ID = "0";

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

void BinStats::WriteToBinningFile(const std::string& binning_file, const ScaffoldsPaths &scaffolds_paths) {
    std::ofstream out(binning_file);

    for (const auto &path_entry : scaffolds_paths) {
      const std::string& scaffold_name = path_entry.first;

      std::vector<EdgeId> scaffold_path(path_entry.second.begin(), path_entry.second.end());
      BinId new_bin_id = ChooseMajorBin(scaffold_path);
      out << scaffold_name << "\t" << (new_bin_id == UNBINNED ? UNBINNED_ID : bin_labels_.at(new_bin_id)) << "\n";
    }
}

BinStats::BinId BinStats::ChooseMajorBin(const std::vector<debruijn_graph::EdgeId>& path) {
    std::vector<size_t> bins_lengths(bins_.size());
    size_t max_len = 0;
    BinId major_bin = UNBINNED;
    for (EdgeId edge : path) {
        if (unbinned_edges_.count(edge) > 0)
            continue;

        size_t length = graph_.length(edge);
        const auto& bins = edges_binning_.at(edge);
        for (auto bin_id : bins) {
            size_t& cur_bin_len = bins_lengths[bin_id];
            cur_bin_len += length;
            if (cur_bin_len > max_len) {
                max_len = cur_bin_len;
                major_bin = bin_id;
            }
        }
    }

    return major_bin;
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
}
