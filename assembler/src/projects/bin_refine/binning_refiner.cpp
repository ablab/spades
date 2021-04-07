//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "binning.hpp"
#include "labels_propagation.hpp"
#include "labels_correction.hpp"
#include "majority_length_strategy.hpp"
#include "max_likelihood_strategy.hpp"

#include "assembly_graph/core/graph.hpp"
#include "pipeline/graph_pack.hpp"
#include "toolchain/utils.hpp"
#include "utils/segfault_handler.hpp"

#include <clipp/clipp.h>

using namespace debruijn_graph;
using namespace bin_stats;

enum class AssignStrategy {
    MajorityLength,
    MaxLikelihood
};

// TODO: add opportunity to choose few types to run in some order
enum class RefinerType {
    Propagation,
    Correction
};

struct gcfg {
    size_t k = 55;
    std::string graph;
    std::string binning_file;
    std::string output_file;
    AssignStrategy assignment_strategy = AssignStrategy::MajorityLength;
    double eps = 0.01;
    double labeled_alpha = 0.6;
    double unlabeled_alpha = 0.999;
    bool allow_multiple = false;
    RefinerType refiner_type = RefinerType::Propagation;
};

static void process_cmdline(int argc, char** argv, gcfg& cfg) {
  using namespace clipp;

  auto cli = (
      cfg.graph << value("graph (in binary or GFA)"),
      cfg.binning_file << value("file with binning from binner in .tsv format"),
      cfg.output_file << value("path to file to write binning after propagation"),
      (option("-e") & value("eps", cfg.eps)) % "convergence relative tolerance threshold",
      (option("-m").set(cfg.allow_multiple) % "allow multiple bin assignment"),
      (with_prefix("-S",
                   option("max").set(cfg.assignment_strategy, AssignStrategy::MajorityLength) |
                   option("mle").set(cfg.assignment_strategy, AssignStrategy::MaxLikelihood)) % "binning assignment strategy"),
      in_sequence(option("-Rcorr").set(cfg.refiner_type, RefinerType::Correction),
           (option("-la") & value("labeled alpha", cfg.labeled_alpha)) % "labels correction alpha for labeled data",
           (option("-ua") & value("unlabeled alpha", cfg.unlabeled_alpha)) % "labels correction alpha for unlabeled data") |
      option("-Rprop").set(cfg.refiner_type, RefinerType::Propagation) % "binning refiner type"
  );

  auto result = parse(argc, argv, cli);
  if (!result) {
    std::cout << make_man_page(cli, argv[0]);
    exit(1);
  }
}

std::unique_ptr<BinningAssignmentStrategy> get_strategy(const gcfg &cfg) {
    switch (cfg.assignment_strategy) {
        default:
            FATAL_ERROR("Unknown binning assignment strategy");
        case AssignStrategy::MajorityLength:
            return std::make_unique<MajorityLengthBinningAssignmentStrategy>(cfg.allow_multiple);
        case AssignStrategy::MaxLikelihood:
            return std::make_unique<MaxLikelihoodBinningAssignmentStrategy>(cfg.allow_multiple);
    }
}

std::unique_ptr<BinningRefiner> get_refiner(const gcfg& cfg, const Graph& graph, size_t num_bins) {
    switch (cfg.refiner_type) {
        default:
            FATAL_ERROR("Unknown binning refiner type");
        case RefinerType::Propagation:
            return std::make_unique<LabelsPropagation>(graph, num_bins, cfg.eps);
        case RefinerType::Correction:
            return std::make_unique<LabelsCorrection>(graph, num_bins, cfg.eps, cfg.labeled_alpha, cfg.unlabeled_alpha);
    }
}

int main(int argc, char** argv) {
  utils::segfault_handler sh;
  gcfg cfg;

  process_cmdline(argc, argv, cfg);

  toolchain::create_console_logger();

  START_BANNER("Binning refiner & propagator");

  try {
      auto assignment_strategy = get_strategy(cfg);

      std::unique_ptr<io::IdMapper<std::string>> id_mapper(new io::IdMapper<std::string>());

      gfa::GFAReader gfa(cfg.graph);
      INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links() << ", paths: " << gfa.num_paths());
      VERIFY_MSG(gfa.k() != -1U, "Failed to determine k-mer length");
      VERIFY_MSG(gfa.k() % 2 == 1, "k-mer length must be odd");

      debruijn_graph::GraphPack gp(gfa.k(), "", 0);
      gfa.to_graph(gp.get_mutable<Graph>(), id_mapper.get());

      const auto& graph = gp.get<Graph>();
      INFO("Graph loaded. Total vertices: " << graph.size() << ", total edges: " << graph.e_size());

      // TODO: For now the edges is a set, we need to decide what to do with
      // repeats (so, we may want to count multiplicity here somehow)
      BinStats::ScaffoldsPaths scaffolds_paths;
      for (const auto &path : gfa.paths()) {
          const std::string &name = path.name;
          // SPAdes outputs paths of scaffolds in the file, so we need to strip the path segment id from the end
          scaffolds_paths[name.substr(0, name.find_last_of('_'))].insert(path.edges.begin(), path.edges.end());
      }

      BinStats binning(graph);
      binning.LoadBinning(cfg.binning_file, scaffolds_paths);

      INFO("Initial binning:\n" << binning);
      auto binning_refiner = get_refiner(cfg, graph, binning.bins().size());
      auto soft_edge_labels = binning_refiner->RefineBinning(binning);
      INFO("Assigning edges to bins");
      binning.AssignEdgeBins(soft_edge_labels, *assignment_strategy);
      INFO("Final binning:\n" << binning);
      INFO("Writing final binning");
      binning.WriteToBinningFile(cfg.output_file, scaffolds_paths,
                                 soft_edge_labels, *assignment_strategy,
                                 *id_mapper);
  } catch (const std::string& s) {
      std::cerr << s << std::endl;
      return EINTR;
  } catch (const std::exception& e) {
      std::cerr << "ERROR: " << e.what() << std::endl;
      return EINTR;
  }
  INFO("Binning refining & propagation finished. Thanks for useful refining!");

  return 0;
}
