//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"
#include "binning.hpp"
#include "binning_propagation.hpp"

#include "modules/alignment/kmer_mapper.hpp"
#include "pipeline/graph_pack.hpp"
#include "toolchain/utils.hpp"
#include "utils/segfault_handler.hpp"

#include <clipp/clipp.h>
#include <iostream>
#include <string>

using namespace debruijn_graph;
using namespace bin_stats;

struct gcfg {
  size_t k = 55;
  std::string graph;
  std::string binning_file;
  std::string scaffolds_file;
  std::string output_file;
  double eps = 0.01;
};

static void process_cmdline(int argc, char** argv, gcfg& cfg) {
  using namespace clipp;

  auto cli = (
      cfg.graph << value("graph (in binary or GFA)"),
      cfg.binning_file << value("file with binning from binner in .tsv format"),
      cfg.scaffolds_file << value("scaffolds in .fasta format"),
      cfg.output_file << value("path to file to write binning after propagation"),
      (option("-e") & value("eps", cfg.eps)) % "iteration min epsilon"
  );

  auto result = parse(argc, argv, cli);
  if (!result) {
    std::cout << make_man_page(cli, argv[0]);
    exit(1);
  }
}

std::vector<EdgeId> conjugate_path(const std::vector<EdgeId> &path,
                                   const debruijn_graph::ConjugateDeBruijnGraph &g) {
    std::vector<EdgeId> cpath;
    for (auto it = path.crbegin(), e = path.crend(); it != e; ++it) {
        cpath.push_back(g.conjugate(*it));
    }
    return cpath;
}

int main(int argc, char** argv) {
  utils::segfault_handler sh;
  gcfg cfg;

  process_cmdline(argc, argv, cfg);

  toolchain::create_console_logger();

  START_BANNER("Binning refiner & propagator");

  try {
      std::unique_ptr<io::IdMapper<std::string>> id_mapper(new io::IdMapper<std::string>());

      gfa::GFAReader gfa(cfg.graph);
      INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links());
      VERIFY_MSG(gfa.k() != -1U, "Failed to determine k-mer length");
      VERIFY_MSG(gfa.k() % 2 == 1, "k-mer length must be odd");

      debruijn_graph::GraphPack gp(gfa.k(), "", 0);
      gfa.to_graph(gp.get_mutable<Graph>(), id_mapper.get());

      const auto& graph = gp.get<Graph>();
      INFO("Graph loaded. Total vertices: " << graph.size() << " Total edges: " << graph.e_size());

      // TODO: For now the edges is a set, we need to decide what to do with
      // repeats (so, we may want to count multiplicity here somehow)
      BinStats::ScaffoldsPaths scaffolds_paths;
      for (const auto &path : gfa.paths()) {
          const std::string &name = path.name;
          // SPAdes outputs paths of scaffolds in the file, so we need to strip the path segment id from the end
          scaffolds_paths[name.substr(0, name.find_last_of("_"))].insert(path.edges.begin(), path.edges.end());
      }

      BinStats binning(graph);
      binning.LoadBinning(cfg.binning_file, cfg.scaffolds_file,
                          scaffolds_paths);

      INFO("" << binning);
      BinningPropagation(graph, cfg.eps).PropagateBinning(binning);
      INFO("" << binning);
      binning.WriteToBinningFile(cfg.output_file, cfg.scaffolds_file, scaffolds_paths);

  } catch (const std::string& s) {
      std::cerr << s << std::endl;
      return EINTR;
  } catch (const std::exception& e) {
      std::cerr << "ERROR: " << e.what() << std::endl;
      return EINTR;
  }

  return 0;
}
