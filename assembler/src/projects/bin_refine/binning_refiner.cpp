//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "alpha_propagation.hpp"
#include "binning.hpp"
#include "labels_propagation.hpp"
#include "link_index.hpp"
#include "majority_length_strategy.hpp"
#include "max_likelihood_strategy.hpp"
#include "paired_end.hpp"

#include "assembly_graph/core/graph.hpp"
#include "pipeline/config_struct.hpp"
#include "toolchain/utils.hpp"
#include "io/graph/gfa_reader.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "utils/parallel/openmp_wrapper.h"
#include "utils/segfault_handler.hpp"
#include "utils/stl_utils.hpp"
#include "utils/verify.hpp"

#include <clipp/clipp.h>
#include <cstdint>
#include <string>

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
    std::string path_file;
    unsigned nthreads = (omp_get_max_threads() / 2 + 1);
    std::string file = "";
    std::string tmpdir = "tmp";
    unsigned libindex = -1u;
    AssignStrategy assignment_strategy = AssignStrategy::MajorityLength;
    double eps = 1e-5;
    double labeled_alpha = 0.6;
    bool allow_multiple = false;
    RefinerType refiner_type = RefinerType::Propagation;
    bool bin_load = false;
    bool debug = false;
    bool bin_dist = false;
    uint64_t out_options = 0;
};

static void process_cmdline(int argc, char** argv, gcfg& cfg) {
  using namespace clipp;

  auto cli = (
      cfg.graph << value("graph (in binary or GFA)"),
      cfg.binning_file << value("file with binning from binner in .tsv format"),
      cfg.output_file << value("path to file to write binning after propagation"),
      (option("--paths") & value("contig.paths", cfg.path_file)) % "use contig paths from file",
      (option("--dataset") & value("yaml", cfg.file)) % "dataset description (in YAML)",
      (option("-l") & integer("value", cfg.libindex)) % "library index (0-based, default: 0)",
      (option("-t") & integer("value", cfg.nthreads)) % "# of threads to use",
      (option("--tmp-dir") & value("dir", cfg.tmpdir)) % "scratch directory to use",
      (option("-e") & value("eps", cfg.eps)) % "convergence relative tolerance threshold",
      (option("-m").set(cfg.allow_multiple) % "allow multiple bin assignment"),
      (with_prefix("-S",
                   option("max").set(cfg.assignment_strategy, AssignStrategy::MajorityLength) |
                   option("mle").set(cfg.assignment_strategy, AssignStrategy::MaxLikelihood)) % "binning assignment strategy"),
      (with_prefix("-R",
                   option("corr").set(cfg.refiner_type, RefinerType::Correction) |
                   option("prop").set(cfg.refiner_type, RefinerType::Propagation)) % "binning refiner type"),
      (option("--cami").call([&]{ cfg.out_options |= OutputOptions::CAMI; }) % "use CAMI bioboxes binning format"),
      (option("--zero-bin").call([&] { cfg.out_options |= OutputOptions::EmitZeroBin; }) % "emit zero bin for unbinned sequences"),
      (option("--tall-multi").call([&] { cfg.out_options |= OutputOptions::TallMulti; }) % "use tall table for multiple binning result"),
      (option("--bin-dist").set(cfg.bin_dist) % "estimate pairwise bin distance (could be slow on large graphs!)"),
      (option("-la") & value("labeled alpha", cfg.labeled_alpha)) % "labels correction alpha for labeled data",
      (option("--bin-load").set(cfg.bin_load)) % "load binary-converted reads from tmpdir (developer option)",
      (option("--debug").set(cfg.debug)) % "produce lots of debug data (developer option)"
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

std::unique_ptr<AlphaAssigner> get_alpha_assigner(const gcfg &cfg,
                                                  const binning::LinkIndex &links,
                                                  const Graph &graph) {
    switch (cfg.refiner_type) {
        default:
            FATAL_ERROR("Unknown binning refiner type");
        case RefinerType::Propagation:
            return std::make_unique<PropagationAssigner>(graph);
        case RefinerType::Correction:
            //todo different choice for metaalpha?
            return std::make_unique<AlphaPropagator>(graph, links, cfg.labeled_alpha);
    }
}

//std::unique_ptr<BinningRefiner> get_refiner(const gcfg &cfg,
//                                            const binning::LinkIndex &links,
//                                            const Graph &graph) {
//    switch (cfg.refiner_type) {
//        default:
//            FATAL_ERROR("Unknown binning refiner type");
//        case RefinerType::Propagation:
//            return std::make_unique<LabelsPropagation>(graph, links, cfg.eps);
//        case RefinerType::Correction:
//            return std::make_unique<LabelsPropagation>(graph, links, cfg.eps, cfg.labeled_alpha);
//    }
//}

int main(int argc, char** argv) {
  utils::segfault_handler sh;
  gcfg cfg;

  srand(42);
  srandom(42);

  process_cmdline(argc, argv, cfg);

  toolchain::create_console_logger();

  START_BANNER("Binning refiner & propagator");

  cfg.nthreads = spades_set_omp_threads(cfg.nthreads);
  INFO("Maximum # of threads to use (adjusted due to OMP capabilities): " << cfg.nthreads);

  fs::CheckFileExistenceFATAL(cfg.graph);
  fs::CheckFileExistenceFATAL(cfg.binning_file);

  fs::make_dir(cfg.tmpdir);

  try {
      auto assignment_strategy = get_strategy(cfg);

      std::unique_ptr<io::IdMapper<std::string>> id_mapper(new io::IdMapper<std::string>());

      DataSet dataset;
      if (cfg.file != "") {
        dataset.load(cfg.file);

        if (cfg.libindex == -1u)
            cfg.libindex = 0;
        CHECK_FATAL_ERROR(cfg.libindex < dataset.lib_count(), "invalid library index");
      }

      gfa::GFAReader gfa(cfg.graph);
      INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links() << ", paths: " << gfa.num_paths());
      VERIFY_MSG(gfa.k() != -1U, "Failed to determine k-mer length");
      VERIFY_MSG(gfa.k() % 2 == 1, "k-mer length must be odd");

      debruijn_graph::Graph graph(gfa.k());
      gfa.to_graph(graph, id_mapper.get());
      INFO("Graph loaded. Total vertices: " << graph.size() << ", total edges: " << graph.e_size());

      INFO("Gathering edge links");
      binning::GraphLinkIndex links(graph);

      // TODO: For now the edges is a set, we need to decide what to do with
      // repeats (so, we may want to count multiplicity here somehow)
      INFO("Processing scaffolds")
      Binning binning(graph);
      {
          std::vector<Binning::Scaffold> scaffolds_paths;
          size_t slinks = 0;

          if (fs::FileExists(cfg.path_file)) {
              INFO("Using contig paths from " << cfg.path_file);
              std::ifstream pathfile_reader(cfg.path_file);
              const std::string numeric = "0123456789";
              for (std::string scaffold_name; std::getline(pathfile_reader, scaffold_name, '\n'); ) {
                  std::string segment;
                  DEBUG("Processing " << scaffold_name);
                  if (!utils::ends_with(scaffold_name, "'")) {
                      EdgeId last;
                      scaffolds_paths.emplace_back(scaffold_name, Binning::ScaffoldPath{});
                      do {
                          std::getline(pathfile_reader, segment, '\n');
                          DEBUG("Processing segment: " << segment);
                          std::vector<EdgeId> edges;
                          utils::split(segment, ',', std::back_inserter(edges),
                                       [&graph, &numeric, &id_mapper](const std::string &edgestr) {
                                           std::size_t idx = edgestr.find_first_not_of(numeric);
                                           VERIFY(idx != std::string::npos);

                                           EdgeId e = (*id_mapper)[edgestr.substr(0, idx)];
                                           if (edgestr[idx] == '-')
                                               e = graph.conjugate(e);
                                           else {
                                               VERIFY(edgestr[idx] == '+');
                                           }

                                           return e;
                                       });

                          if (last != EdgeId()) {
                              // If this is a proper scaffold (multiple paths), then link
                              // paths as there might be no graph connectivity
                              links.add(last, edges.front());
                              slinks += 1;
                          }
                          last = edges.back();
                          scaffolds_paths.back().second.insert(edges.begin(), edges.end());
                      } while (utils::ends_with(segment, ";"));
                  }  else {
                      // Skip conjugate path
                      DEBUG("Skipping entries");
                      do {
                          std::getline(pathfile_reader, segment, '\n');
                      } while (utils::ends_with(segment, ";"));
                  }
              }
          } else {
              INFO("Using scaffold paths from GFA");

              std::string scaffold_name;
              EdgeId last;
              for (const auto &path : gfa.paths()) {
                  const std::string &name = path.name;
                  std::string cname = name.substr(0, name.find_last_of('_'));
                  // SPAdes outputs paths of scaffolds in the file, so we need to strip the path segment id from the end
                  if (cname != scaffold_name) {
                      scaffold_name = cname;
                      scaffolds_paths.emplace_back(scaffold_name, Binning::ScaffoldPath{});
                  } else if (last != EdgeId()) {
                      // If this is a proper scaffold (multiple paths), then link
                      // paths as there might be no graph connectivity
                      links.add(last, path.edges.front());
                      slinks += 1;
                  }
                  last = path.edges.back();
                  scaffolds_paths.back().second.insert(path.edges.begin(), path.edges.end());
              }
          }

          INFO("Added " << slinks << " additional scaffold links");
          binning.InitScaffolds(scaffolds_paths);
          INFO("Processed total " << scaffolds_paths.size() << " scaffold paths");
      }

      binning::LinkIndex pe_links(graph);
      if (cfg.libindex != -1u) {
          INFO("Processing paired-end reads");
          debruijn_graph::config::init_libs(dataset, cfg.nthreads, cfg.tmpdir);

          auto &lib = dataset[cfg.libindex];
          if (lib.is_paired()) {
              binning::FillPairedEndLinks(pe_links, lib, graph,
                                          cfg.tmpdir, cfg.nthreads, cfg.bin_load, cfg.debug);

              for (auto it = pe_links.begin(), end = pe_links.end(); it != end; ++it) {
                  EdgeId e1 = it.key();

                  for (const auto &link: it.value()) {
                      if (link.w == 1)
                          continue;
                      if (!(e1 <= link.e)) // do not process the same link twice
                          continue;

                      links.increment(e1, link.e, log2(link.w));
                  }
              }
          } else {
              WARN("Only paired-end libraries are supported for links");
          }
      }

      binning.LoadBinning(cfg.binning_file, cfg.out_options & OutputOptions::CAMI);
      INFO("Initial binning:\n" << binning);

      auto alpha_assigner = get_alpha_assigner(cfg, links, graph);
      auto alpha_assignment = alpha_assigner->GetAlphaAssignment(binning);
      auto binning_refiner = std::make_unique<LabelsPropagation>(graph, links, alpha_assignment, cfg.eps);
//      auto binning_refiner = get_refiner(cfg, links, graph);
      auto soft_edge_labels = binning_refiner->RefineBinning(binning);

      INFO("Assigning edges & scaffolds to bins");
      binning.AssignBins(soft_edge_labels, *assignment_strategy);
      INFO("Final binning:\n" << binning);
      INFO("Estimating pairwise bin distances");
      auto dist = binning.BinDistance(soft_edge_labels);
      {
          size_t nbins = binning.bins().size();

          std::ofstream out_dist(cfg.output_file + ".bin_dist");
          for (size_t i = 0; i < nbins; ++i)
              out_dist << '\t' << binning.bin_labels().at(i);
          out_dist << std::endl;

          for (size_t i = 0; i < nbins; ++i) {
              out_dist << binning.bin_labels().at(i);
              for (size_t j = 0; j < nbins; ++j)
                  out_dist << '\t' << dist(i, j);
              out_dist << std::endl;
          }
      }
      INFO("Writing final binning");
      binning.WriteToBinningFile(cfg.output_file, cfg.out_options,
                                 soft_edge_labels, *assignment_strategy,
                                 *id_mapper);
      if (cfg.debug) {
          INFO("Dumping links");
          pe_links.dump(cfg.output_file + ".pe_links", *id_mapper);
          links.dump(cfg.output_file + ".graph_links", *id_mapper);
      }

      if (!cfg.debug)
          fs::remove_dir(cfg.tmpdir);
  } catch (const std::string& s) {
      std::cerr << s << std::endl;
      return EINTR;
  } catch (const std::exception& e) {
      std::cerr << "ERROR: " << e.what() << std::endl;
      return EINTR;
  }
  INFO("Binning refining & propagation finished. Thanks for interesting data!");

  return 0;
}
