//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2021-2022 Saint Petersburg State University
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
#include "read_splitting.hpp"

#include "assembly_graph/core/graph.hpp"
#include "configs/config_struct.hpp"
#include "toolchain/utils.hpp"
#include "io/graph/gfa_reader.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "utils/logger/logger.hpp"
#include "utils/parallel/openmp_wrapper.h"
#include "utils/segfault_handler.hpp"
#include "utils/stl_utils.hpp"
#include "utils/verify.hpp"

#include <clipp/clipp.h>
#include <cstdint>
#include <filesystem>

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
    std::filesystem::path graph;
    std::filesystem::path binning_file;
    std::filesystem::path prefix;
    std::filesystem::path path_file;
    unsigned nthreads = (omp_get_max_threads() / 2 + 1);
    std::string file = "";
    std::filesystem::path tmpdir = "tmp";
    unsigned libindex = -1u;
    AssignStrategy assignment_strategy = AssignStrategy::MaxLikelihood;
    double eps = 1e-5;
    unsigned niter = 5000;
    double labeled_alpha = 0.6;
    bool no_unbinned_bin = false;
    bool alpha_propagation = false;
    double metaalpha = 0.6;
    size_t length_threshold = 5000;
    size_t distance_bound = 10000;
    bool allow_multiple = false;
    RefinerType refiner_type = RefinerType::Correction;
    bool bin_load = false;
    bool debug = false;
    bool bin_dist = false;
    uint64_t out_options = 0;
    bool split_reads = false;
    double bin_weight_threshold = 0.1;
};

static void process_cmdline(int argc, char** argv, gcfg& cfg) {
  using namespace clipp;

  std::string path_file, graph, binning_file, prefix, tmpdir = "tmp";
  auto cli = (
      graph << value("graph (in binary or GFA)"),
      binning_file << value("file with binning from binner in .tsv format"),
      prefix << value("output path to write binning results after propagation"),
      (option("--paths") & value("contig.paths", path_file)) % "use contig paths from file",
      (option("--dataset") & value("yaml", cfg.file)) % "dataset description (in YAML)",
      (option("-l") & integer("value", cfg.libindex)) % "library index (0-based, default: 0)",
      (option("-t") & integer("value", cfg.nthreads)) % "# of threads to use",
      (option("-e") & value("eps", cfg.eps)) % "convergence relative tolerance threshold",
      (option("-n") & integer("value", cfg.niter)) % "maximum number of iterations",
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
      "Sparse propagation options:" % (
          (option("--sparse-propagation").set(cfg.alpha_propagation)) % "Gradually reduce alpha from binned to unbinned edges",
          (option("--no-unbinned-bin").set(cfg.no_unbinned_bin)) % "Do not create a special bin for unbinned contigs",
          (option("-ma") & value("--metaalpha", cfg.metaalpha)) % "Labels correction alpha for sparse propagation procedure",
          (option("-lt") & value("--length-threshold", cfg.length_threshold)) % "Binning will not be propagated to edges longer than threshold",
          (option("-db") & value("--distance-bound", cfg.distance_bound)) % "Binning will not be propagated further than bound"
      ),
      "Read splitting options:" % (
          (option("-r", "--reads").set(cfg.split_reads) % "split reads according to binning"),
          (option("-b", "--bin-weight") & value("threshold", cfg.bin_weight_threshold)) % "reads bin weight threshold"
      ),
      "Developer options:" % (
          (option("--bin-load").set(cfg.bin_load)) % "load binary-converted reads from tmpdir",
          (option("--debug").set(cfg.debug)) % "produce lots of debug data",
          (option("--tmp-dir") & value("dir", tmpdir)) % "scratch directory to use"
      )
  );

  auto result = parse(argc, argv, cli);
  if (!result) {
      std::cout << make_man_page(cli, argv[0]);
      exit(1);
  }

  cfg.prefix = prefix;
  cfg.graph = graph;
  cfg.path_file = path_file;
  cfg.binning_file = binning_file;
  cfg.tmpdir = tmpdir;
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
                                                  const Graph &graph,
                                                  const bin_stats::Binning &binning) {
    switch (cfg.refiner_type) {
        default:
            FATAL_ERROR("Unknown binning refiner type");
        case RefinerType::Propagation:
            return std::make_unique<PropagationAssigner>(graph);
        case RefinerType::Correction:
            if (not cfg.alpha_propagation)
                return std::make_unique<CorrectionAssigner>(graph, cfg.labeled_alpha);

            auto alpha_mask = AlphaPropagator(graph, links, cfg.metaalpha, cfg.eps, cfg.niter,
                                              cfg.length_threshold, cfg.distance_bound,
                                              cfg.prefix / "alpha_stats.tsv").GetAlphaMask(binning);
            return std::make_unique<bin_stats::CorrectionAssigner>(graph, alpha_mask, cfg.labeled_alpha);
    }
}

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

  if (!std::filesystem::exists(cfg.graph))
      FATAL_ERROR("Graph file: " << cfg.graph << " does not exist");

  if (!std::filesystem::exists(cfg.binning_file))
      FATAL_ERROR("Input binning file " << cfg.binning_file << " does not exist");

  if (!std::filesystem::create_directories(cfg.tmpdir))
      FATAL_ERROR("Failed to create temporary directory: " << cfg.tmpdir);

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

      debruijn_graph::Graph graph(0);

      gfa::GFAReader gfa(cfg.graph);
      unsigned gfa_k = gfa.to_graph(graph, id_mapper.get());

      INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links() <<
           ", paths: " << gfa.num_paths() << ", jumps: " << gfa.num_gaplinks());
      if (gfa_k != -1U) {
          INFO("Detected k: " << gfa_k);
          VERIFY_MSG(gfa_k == 0 || gfa_k % 2 == 1, "k-mer length must be odd");
          VERIFY(graph.k() == gfa_k);
      } else {
          INFO("Graph seems to be multiplexed");
      }

      INFO("Graph loaded. Total vertices: " << graph.size() << ", total edges: " << graph.e_size());

      INFO("Gathering edge links");
      binning::GraphLinkIndex links(graph);

      INFO("Adding jump links");
      for (const auto &jump : gfa.jumps())
          links.add(jump.first, jump.second);

      // TODO: For now the edges is a set, we need to decide what to do with
      // repeats (so, we may want to count multiplicity here somehow)
      INFO("Processing scaffolds")
      Binning binning(graph);
      {
          std::vector<Binning::Scaffold> scaffolds_paths;
          size_t slinks = 0;

          if (std::filesystem::exists(cfg.path_file)) {
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
                          utils::split(segment, ",", std::back_inserter(edges),
                                       [&graph, &numeric, &id_mapper](const std::string_view edgestr) {
                                           std::size_t idx = edgestr.find_first_not_of(numeric);
                                           VERIFY(idx != std::string::npos);

                                           EdgeId e = (*id_mapper)[std::string(edgestr.substr(0, idx))];
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
                  // SPAdes outputs paths of scaffolds in the file, so we need
                  // to strip the path segment id from the end However, this
                  // only applies to scaffold names emitted in GFA v1.1 graphs
                  // and not for GFA v1.2
                  std::string cname;
                  if (utils::starts_with(name, "NODE_")) {
                      auto upos = name.find_last_of('_');
                      // GFA v1.1 format: NODE_123_length_56_cov_2.000000_0
                      // GFA v1.2 format: NODE_123_length_56_cov_2.000000

                      // As coverage is always float we also checking for
                      // position of dot symbol. For GFA v1.2 names it will be
                      // after the last underscore and before for GFA v1.1.
                      auto dpos = name.find_last_of('.');
                      cname = dpos < upos ? name.substr(0, upos) : name;
                  } else
                      cname = name;
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
      if (cfg.libindex != -1u || cfg.split_reads)
          debruijn_graph::config::init_libs(dataset, cfg.nthreads, cfg.tmpdir);

      if (cfg.libindex != -1u) {
          INFO("Processing paired-end reads");
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

      bool add_unbinned_bin = !cfg.no_unbinned_bin && cfg.alpha_propagation;
      binning.LoadBinning(cfg.binning_file, cfg.out_options & OutputOptions::CAMI, add_unbinned_bin);
      INFO("Initial binning:\n" << binning);

      auto alpha_assigner = get_alpha_assigner(cfg, links, graph, binning);
      LabelInitializer label_initializer(graph, cfg.length_threshold, add_unbinned_bin);
      auto origin_state = label_initializer.InitLabels(binning);
      auto alpha_assignment = alpha_assigner->GetAlphaAssignment(origin_state);

      const size_t unbinned_length_threshold = cfg.length_threshold;
      std::unordered_set<EdgeId> nonpropagating_edges;
      if (add_unbinned_bin) {
          for (const EdgeId &edge: binning.unbinned_edges()) {
              if (graph.length(edge) >= unbinned_length_threshold) {
                  nonpropagating_edges.insert(edge);
              }
          }
      }
      auto binning_refiner = std::make_unique<LabelsPropagation>(graph, links, alpha_assignment, nonpropagating_edges,
                                                                 cfg.eps, cfg.niter);
      auto soft_edge_labels = binning_refiner->RefineBinning(origin_state);

      INFO("Assigning edges & scaffolds to bins");
      binning.AssignBins(soft_edge_labels, *assignment_strategy);
      INFO("Final binning:\n" << binning);

      if (!std::filesystem::create_directories(cfg.prefix))
          FATAL_ERROR("Failed to create output dir: " << cfg.prefix);

      INFO("Writing final binning");
      binning.WriteToBinningFile(cfg.prefix, cfg.out_options,
                                 soft_edge_labels, *assignment_strategy,
                                 *id_mapper);

      if (cfg.bin_dist) {
          INFO("Estimating pairwise bin distances");
          auto dist = binning.BinDistance(soft_edge_labels);
          size_t nbins = binning.bins().size();

          std::ofstream out_dist(cfg.prefix / "bin_dist.tsv");
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

      if (cfg.debug) {
          INFO("Dumping links");
          pe_links.dump(cfg.prefix / "pe_links.tsv", *id_mapper);
          links.dump(cfg.prefix / "graph_links.tsv", *id_mapper);
      }

      if (cfg.split_reads) {
          auto &lib = dataset[0];

          INFO("Splitting reads started");
          binning::SplitAndWriteReads(graph,
                                      lib,
                                      binning,
                                      soft_edge_labels,
                                      *assignment_strategy,
                                      cfg.tmpdir,
                                      cfg.prefix,
                                      cfg.nthreads,
                                      cfg.bin_weight_threshold);
          INFO("Splitting reads ended");
      }

      if (!cfg.debug)
          if (!std::filesystem::remove_all(cfg.tmpdir))
              FATAL_ERROR("Failed to empty temporary directory: " << cfg.tmpdir);
  } catch (const std::string& s) {
      std::cerr << s << std::endl;
      return EINTR;
  } catch (const std::exception& e) {
      std::cerr << "ERROR: " << e.what() << std::endl;
      return EINTR;
  }
  INFO("Binning refining & propagation finished. Let's analyze these MAGs!");

  return 0;
}
