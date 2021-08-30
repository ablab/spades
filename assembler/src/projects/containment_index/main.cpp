//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "barcode_index_construction.hpp"
#include "scaffold_graph_helper.hpp"

#include "io/binary/read_cloud.hpp"
#include "toolchain/utils.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "utils/parallel/openmp_wrapper.h"
#include "utils/segfault_handler.hpp"
#include "utils/verify.hpp"

#include <clipp/clipp.h>

using namespace debruijn_graph;
using namespace cont_index;

struct gcfg {
  size_t k = 55;
  std::string graph;
  std::string output_file;
  unsigned nthreads = (omp_get_max_threads() / 2 + 1);
  std::string file = "";
  std::string tmpdir = "saves";
  unsigned libindex = -1u;
  bool bin_load = false;
  bool debug = false;
};

static void process_cmdline(int argc, char** argv, gcfg& cfg) {
    using namespace clipp;

    auto cli = (
        cfg.graph << value("graph (in binary or GFA)"),
            cfg.output_file << value("path to file to write binning after propagation"),
            (option("--dataset") & value("yaml", cfg.file)) % "dataset description (in YAML)",
            (option("-l") & integer("value", cfg.libindex)) % "library index (0-based, default: 0)",
            (option("-t") & integer("value", cfg.nthreads)) % "# of threads to use",
            (option("--tmp-dir") & value("tmp", cfg.tmpdir)) % "scratch directory to use",
            (option("--bin-load").set(cfg.bin_load)) % "load binary-converted reads from tmpdir (developer option)",
            (option("--debug").set(cfg.debug)) % "produce lots of debug data (developer option)"
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

    srand(42);
    srandom(42);

    process_cmdline(argc, argv, cfg);

    toolchain::create_console_logger();

    START_BANNER("Containment index builder");

    cfg.nthreads = spades_set_omp_threads(cfg.nthreads);
    INFO("Maximum # of threads to use (adjusted due to OMP capabilities): " << cfg.nthreads);

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

        debruijn_graph::config::init_libs(dataset, cfg.nthreads, cfg.tmpdir);
        barcode_index::FrameBarcodeIndex<debruijn_graph::Graph> barcode_index(graph, 2000);

        auto &lib = dataset[cfg.libindex];
        if (lib.type() == io::LibraryType::Clouds10x) {
            cont_index::ConstructBarcodeIndex(barcode_index, lib, graph, cfg.tmpdir, cfg.nthreads, cfg.bin_load, cfg.debug);
        } else {
            WARN("Only read cloud libraries with barcode tags are supported for links");
        }

        INFO(barcode_index.GetNumberOfBarcodes() << " barcodes in index");
        using BarcodeExtractor = barcode_index::FrameBarcodeIndexInfoExtractor;
        auto barcode_extractor_ptr = std::make_shared<BarcodeExtractor>(barcode_index, graph);
        LinkIndexGraphConstructor link_index_constructor(graph, barcode_extractor_ptr, cfg.nthreads);
        auto scaffold_graph = link_index_constructor.ConstructGraph();

        std::ofstream os(cfg.output_file);
        os << "FirstId\tSecondId\tWeight\n";
        for (const scaffold_graph::ScaffoldGraph::ScaffoldEdge &edge: scaffold_graph.edges()) {
            os << (*id_mapper)[edge.getStart().int_id()] << "\t" << (*id_mapper)[edge.getEnd().int_id()] << "\t" << edge.getWeight() << "\n";
        }

    }

}