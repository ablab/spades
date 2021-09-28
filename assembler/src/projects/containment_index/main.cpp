//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "barcode_index_construction.hpp"
#include "scaffold_graph_helper.hpp"
#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"

#include "io/binary/read_cloud.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/edge_cluster_extractor.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/initial_cluster_storage_builder.hpp"
#include "modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_polisher.hpp"
#include "modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/path_cluster_helper.hpp"
#include "toolchain/utils.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "utils/parallel/openmp_wrapper.h"
#include "utils/segfault_handler.hpp"
#include "utils/verify.hpp"

#include <clipp/clipp.h>

using namespace debruijn_graph;
using namespace cont_index;
using namespace path_extend::read_cloud;

struct gcfg {
  size_t k = 55;
  std::string graph;
  std::string output_dir;
  unsigned nthreads = (omp_get_max_threads() / 2 + 1);
  std::string file = "";
  std::string tmpdir = "saves";
  unsigned libindex = -1u;
  double score_threshold = 3.0;
  bool bin_load = false;
  bool debug = false;
};

static void process_cmdline(int argc, char** argv, gcfg& cfg) {
    using namespace clipp;

    auto cli = (
        cfg.graph << value("graph (in binary or GFA)"),
            cfg.output_dir << value("path to output directory"),
            (option("--dataset") & value("yaml", cfg.file)) % "dataset description (in YAML)",
            (option("-l") & integer("value", cfg.libindex)) % "library index (0-based, default: 0)",
            (option("-t") & integer("value", cfg.nthreads)) % "# of threads to use",
            (option("--tmp-dir") & value("tmp", cfg.tmpdir)) % "scratch directory to use",
            (option("--bin-load").set(cfg.bin_load)) % "load binary-converted reads from tmpdir (developer option)",
            (option("--debug").set(cfg.debug)) % "produce lots of debug data (developer option)",
            (option("--score") & value("score", cfg.score_threshold)) % "Score threshold for link index"
    );

    auto result = parse(argc, argv, cli);
    if (!result) {
        std::cout << make_man_page(cli, argv[0]);
        exit(1);
    }
}

struct TimeTracerRAII {
  TimeTracerRAII(llvm::StringRef program_name,
                 unsigned granularity = 500,
                 const std::string &prefix = "", const std::string &suffix = "") {
      time_trace_file_ = prefix + "time_trace_" + suffix + ".json";
      llvm::timeTraceProfilerInitialize(granularity, program_name);
  }
  ~TimeTracerRAII() {
      if (auto E = llvm::timeTraceProfilerWrite(time_trace_file_, "cont-index")) {
          handleAllErrors(std::move(E),
                          [&](const llvm::StringError &SE) {
                            ERROR("" << SE.getMessage() << "\n");
                          });
          return;
      } else {
          INFO("Time trace is written to: " << time_trace_file_);
      }
      llvm::timeTraceProfilerCleanup();
  }

  std::string time_trace_file_;
};

int main(int argc, char** argv) {
    utils::segfault_handler sh;
    gcfg cfg;

    srand(42);
    srandom(42);

    process_cmdline(argc, argv, cfg);

    toolchain::create_console_logger();

    START_BANNER("Containment index builder");

    cfg.nthreads = spades_set_omp_threads(cfg.nthreads);
    INFO("Maximum # of threads to use (adjusted due to OMP capabilities): " << cfg.nthreads);

    fs::make_dir(cfg.output_dir);
    fs::make_dir(cfg.tmpdir);

    INFO("Loading graph");
    std::unique_ptr<io::IdMapper<std::string>> id_mapper(new io::IdMapper<std::string>());

    gfa::GFAReader gfa(cfg.graph);
    INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links() << ", paths: " << gfa.num_paths());
    VERIFY_MSG(gfa.k() != -1U, "Failed to determine k-mer length");
//    VERIFY_MSG(gfa.k() % 2 == 1, "k-mer length must be odd");

    debruijn_graph::Graph graph(gfa.k());
    gfa.to_graph(graph, id_mapper.get());
    INFO("Graph loaded. Total vertices: " << graph.size() << ", total edges: " << graph.e_size());

    INFO("Building barcode index");
    if (cfg.libindex != -1u) {
        INFO("Processing paired-end reads");

        DataSet dataset;
        if (cfg.file != "") {
            dataset.load(cfg.file);

            if (cfg.libindex == -1u)
                cfg.libindex = 0;
            CHECK_FATAL_ERROR(cfg.libindex < dataset.lib_count(), "invalid library index");
        }

        //fixme configs
        const size_t frame_size = 2000;
        const size_t read_linkage_distance = 5000;

        const double graph_score_threshold = 3.99;
        const size_t tail_threshold = 20000;
        const size_t length_threshold = 5000;
        const size_t count_threshold = 3;

        const double relative_score_threshold = 10.0;
        const size_t min_read_threshold = 2;

        debruijn_graph::config::init_libs(dataset, cfg.nthreads, cfg.tmpdir);
        barcode_index::FrameBarcodeIndex<debruijn_graph::Graph> barcode_index(graph, read_linkage_distance);

        std::unique_ptr<TimeTracerRAII> traceraii;
        traceraii.reset(new TimeTracerRAII(argv[0],500));
        INFO("Time tracing is enabled");

        TIME_TRACE_SCOPE("Containment index");

        auto &lib = dataset[cfg.libindex];
        if (lib.type() == io::LibraryType::Clouds10x) {
            cont_index::ConstructBarcodeIndex(barcode_index,
                                              lib,
                                              graph,
                                              cfg.tmpdir,
                                              cfg.nthreads,
                                              frame_size,
                                              cfg.bin_load,
                                              cfg.debug);
        } else {
            WARN("Only read cloud libraries with barcode tags are supported for links");
        }
        INFO("Barcode index size: " << barcode_index.size());
        using BarcodeExtractor = barcode_index::FrameBarcodeIndexInfoExtractor;
        auto barcode_extractor_ptr = std::make_shared<BarcodeExtractor>(barcode_index, graph);

        auto scaffold_graph = cont_index::GetTellSeqScaffoldGraph(graph, barcode_extractor_ptr, graph_score_threshold,
                                                                  length_threshold,
                                                                  tail_threshold, count_threshold, cfg.nthreads,
                                                                  cfg.bin_load,
                                                                  cfg.debug, cfg.output_dir, id_mapper.get());

        INFO("Constructing initial cluster storage");
        path_extend::read_cloud::PathExtractionParams path_extraction_params(read_linkage_distance,
                                                                             relative_score_threshold,
                                                                             min_read_threshold,
                                                                             length_threshold);
        size_t cluster_storage_builder_threads = cfg.nthreads;
        std::set<scaffold_graph::ScaffoldVertex> target_edges;
        std::copy(scaffold_graph.vbegin(), scaffold_graph.vend(), std::inserter(target_edges, target_edges.begin()));
        auto edge_cluster_extractor =
            std::make_shared<cluster_storage::AccurateEdgeClusterExtractor>(graph, barcode_extractor_ptr,
                                                                            read_linkage_distance, min_read_threshold);
        auto storage_builder =
            std::make_shared<cluster_storage::EdgeInitialClusterStorageBuilder>(graph, edge_cluster_extractor,
                                                                                target_edges, read_linkage_distance,
                                                                                min_read_threshold,
                                                                                cluster_storage_builder_threads);
        auto storage =
            std::make_shared<cluster_storage::InitialClusterStorage>(storage_builder->ConstructInitialClusterStorage());
        path_extend::read_cloud::ScaffoldGraphPathClusterHelper path_extractor_helper(graph,
                                                                                      barcode_extractor_ptr,
                                                                                      storage,
                                                                                      read_linkage_distance,
                                                                                      cfg.nthreads);

        INFO(storage->get_cluster_storage().Size() << " initial clusters");
        auto all_clusters = path_extractor_helper.GetAllClusters(scaffold_graph);
        INFO(all_clusters.size() << " total clusters");
        std::map<size_t, size_t> size_hist;
        for (const auto &cluster: all_clusters) {
            size_hist[cluster.Size()]++;
        }
        std::ofstream sizes(fs::append_path(cfg.output_dir, "cluster_sizes.tsv"));
        for (const auto &entry: size_hist) {
            sizes << entry.first << "\t" << entry.second << std::endl;
        }
        auto path_clusters = path_extractor_helper.GetPathClusters(all_clusters);
        INFO(path_clusters.size() << " path clusters");
        size_hist.clear();
        for (const auto &cluster: path_clusters) {
            size_hist[cluster.Size()]++;
        }
        std::ofstream path_sizes(fs::append_path(cfg.output_dir, "path_sizes.tsv"));
        for (const auto &entry: size_hist) {
            path_sizes << entry.first << "\t" << entry.second << std::endl;
        }
    }

}
