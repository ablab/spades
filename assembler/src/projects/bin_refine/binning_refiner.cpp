//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "binning.hpp"

#include "modules/alignment/kmer_mapper.hpp"
#include "pipeline/graph_pack.hpp"
#include "toolchain/utils.hpp"
#include "utils/segfault_handler.hpp"

#include <clipp/clipp.h>
#include <iostream>
#include <string>

using namespace debruijn_graph;

struct gcfg {
  size_t k = 55;
  std::string graph;
  std::string binning_file;
  std::string scaffolds_file;
  double eps = 0.5;
};

static void process_cmdline(int argc, char** argv, gcfg& cfg) {
  using namespace clipp;

  auto cli = (
      cfg.graph << value("graph (in binary or GFA)"),
      cfg.binning_file << value("file with binning from binner in .tsv format"),
      cfg.scaffolds_file << value("scaffolds in .fasta format"),
      (option("-k") & integer("value", cfg.k)) % "k-mer length to use",
      (option("-e") & value("eps", cfg.eps)) % "iteration min epsilon"
  );

  auto result = parse(argc, argv, cli);
  if (!result) {
    std::cout << make_man_page(cli, argv[0]);
    exit(1);
  }
}

int main(int argc, char** argv) {
  utils::segfault_handler sh;
  gcfg cfg;

  process_cmdline(argc, argv, cfg);

  toolchain::create_console_logger();

  try {
      // FIXME: infer k from GFA
      debruijn_graph::GraphPack gp(cfg.k, "", 0);
      const auto& graph = gp.get<Graph>();
      toolchain::LoadGraph(gp, cfg.graph);
      // FIXME: do not need this
      gp.get_mutable<KmerMapper<Graph>>().Attach();
      gp.EnsureBasicMapping();

      bin_stats::BinStats binning(graph);
      binning.LoadBinning(cfg.binning_file, cfg.scaffolds_file, gp);

      INFO("" << binning);
  } catch (const std::string& s) {
      std::cerr << s << std::endl;
      return EINTR;
  } catch (const std::exception& e) {
      std::cerr << "ERROR: " << e.what() << std::endl;
      return EINTR;
  }

  return 0;
}
