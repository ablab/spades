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

const std::string BinStats::SCAFFOLD_NAME_PREFIX = "NODE_";
const std::string BinStats::UNBINNED_ID = "0";


void BinStats::ScaffoldsToEdges(const std::string& scaffolds_file,
                                const debruijn_graph::GraphPack &gp) {
  fs::CheckFileExistenceFATAL(scaffolds_file);

  edges_binning_.clear();
  unbinned_edges_.clear();
  auto mapper = MapperInstance(gp);

  io::SingleStream scaffolds_stream = io::EasyStream(scaffolds_file, false);
  io::SingleRead scaffold;

  while (!scaffolds_stream.eof()) {
    scaffolds_stream >> scaffold;
    const std::string scaffold_name = scaffold.name();
    if (!scaffolds_binning_.count(scaffold_name)) {
      DEBUG("No such scaffold " + scaffold_name + " with name " + scaffold_name);
      continue;
    }
    
    BinId bin_id = scaffolds_binning_[scaffold_name];
    const auto scaffold_path = mapper->MapRead(scaffold);

    for (size_t i = 0; i < scaffold_path.size(); ++i) {
        EdgeId edge = scaffold_path[i].first;
        edges_binning_[edge].insert(bin_id);
        edges_binning_[graph_.conjugate(edge)].insert(bin_id);
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

void BinStats::LoadBinning(const std::string& binning_file, const std::string& scaffolds_file,
                           const debruijn_graph::GraphPack &gp) {
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

  ScaffoldsToEdges(scaffolds_file, gp);
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
