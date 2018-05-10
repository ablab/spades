//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"
#include "utils/extension_index/kmer_extension_index_builder.hpp"

#include "io/reads/read_processor.hpp"
#include "io/reads/io_helper.hpp"
#include "io/graph/gfa_writer.hpp"
#include "io/graph/fastg_writer.hpp"

#include "io/dataset_support/read_converter.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "io/reads/osequencestream.hpp"

#include "assembly_graph/construction/debruijn_graph_constructor.hpp"
#include "common/pipeline/graphio.hpp"

#include "version.hpp"

#include <clipp/clipp.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <string>

using namespace std;
void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

enum class output_type {
    unitigs, fastg, gfa, spades
};

struct gcfg {
    gcfg()
        : k(21), tmpdir("tmp"), outfile("-"),
          nthreads(omp_get_max_threads() / 2 + 1), buff_size(512ULL << 20),
          mode(output_type::unitigs)
    {}

    unsigned k;
    std::string file;
    std::string tmpdir;
    std::string outfile;
    unsigned nthreads;
    size_t buff_size;
    enum output_type mode;
};


void process_cmdline(int argc, char **argv, gcfg &cfg) {
  using namespace clipp;

  auto cli = (
      cfg.file << value("dataset description (in YAML)"),
      cfg.outfile << value("output filename"),
      (option("-k") & integer("value", cfg.k)) % "k-mer length to use",
      (option("-t") & integer("value", cfg.nthreads)) % "# of threads to use",
      (option("-tmpdir") & value("dir", cfg.tmpdir)) % "scratch directory to use",
      (option("-b") & integer("value", cfg.buff_size)) % "sorting buffer size, per thread",
      one_of(option("-unitigs").set(cfg.mode, output_type::unitigs) % "produce unitigs (default)",
             option("-fastg").set(cfg.mode, output_type::fastg) % "produce graph in FASTG format",
             option("-gfa").set(cfg.mode, output_type::gfa) % "produce graph in GFA1 format",
             option("-spades").set(cfg.mode, output_type::spades) % "produce graph in SPAdes internal format")
  );

  auto result = parse(argc, argv, cli);
  if (!result) {
      std::cout << make_man_page(cli, argv[0]);
      exit(1);
  }
}


int main(int argc, char* argv[]) {
    gcfg cfg;

    srand(42);
    srandom(42);

    process_cmdline(argc, argv, cfg);

    try {
        unsigned nthreads = cfg.nthreads;
        unsigned k = cfg.k;
        std::string tmpdir = cfg.tmpdir, dataset_desc = cfg.file;
        size_t buff_size = cfg.buff_size;

        create_console_logger();

        INFO("Starting SPAdes standalone graph builder, built from " SPADES_GIT_REFSPEC ", git revision " SPADES_GIT_SHA1);

        if (k < runtime_k::MIN_K)
            FATAL_ERROR("k-mer size " << k << " is too low");
        if (k >= runtime_k::MAX_K)
            FATAL_ERROR("k-mer size " << k << " is too high, recompile with larger SPADES_MAX_K option");
        if (k % 2 == 0)
            FATAL_ERROR("k-mer size must be odd");


        INFO("K-mer length set to " << k);
        switch (cfg.mode) {
            case output_type::unitigs:
                INFO("Producing unitigs only");
                break;
            case output_type::fastg:
                INFO("Producing graph in FASTG format");
                break;
            case output_type::gfa:
                INFO("Producing graph in GFA1 format");
                break;
            case output_type::spades:
                INFO("Producing graph in SPAdes internal format");
                break;
        }

        nthreads = spades_set_omp_threads(nthreads);
        INFO("Maximum # of threads to use (adjusted due to OMP capabilities): " << nthreads);

        fs::make_dir(tmpdir);
        auto workdir = fs::tmp::make_temp_dir(tmpdir, "construction");

        io::DataSet<debruijn_graph::config::LibraryData> dataset;
        dataset.load(dataset_desc);
        // FIXME: Get rid of this "/" junk
        debruijn_graph::config::init_libs(dataset, nthreads, buff_size, tmpdir + "/");

        std::vector<size_t> libs_for_construction;
        for (size_t i = 0; i < dataset.lib_count(); ++i)
            if (dataset[i].is_graph_contructable())
                libs_for_construction.push_back(i);

        auto read_streams = io::single_binary_readers_for_libs(dataset, libs_for_construction, true, true);

        // Step 1: build extension index
        VERIFY_MSG(read_streams.size(), "No input streams specified");
        utils::DeBruijnExtensionIndex<> ext_index(k);

        utils::DeBruijnExtensionIndexBuilder().BuildExtensionIndexFromStream(workdir,
                                                                             ext_index, read_streams, nullptr, buff_size);

        // Step 2: extract unbranching paths
        bool keep_perfect_loops = true;
        std::vector<Sequence> edge_sequences;
        unsigned nchunks = 16 * omp_get_max_threads();
        if (keep_perfect_loops)
            edge_sequences = debruijn_graph::UnbranchingPathExtractor(ext_index, k).ExtractUnbranchingPathsAndLoops(nchunks);
        else
            edge_sequences = debruijn_graph::UnbranchingPathExtractor(ext_index, k).ExtractUnbranchingPaths(nchunks);

        if (cfg.mode == output_type::unitigs) {
            // Step 3: output stuff

            INFO("Saving unitigs to " << cfg.outfile);
            size_t idx = 1;
            std::ofstream f(cfg.outfile);
            for (const auto &edge: edge_sequences) {
                f << std::string(">") << io::MakeContigId(idx++, edge.size(), "EDGE") << std::endl;
                io::WriteWrapped(edge.str(), f);
            }
        } else {
            // Step 3: build the graph
            INFO("Building graph");
            debruijn_graph::DeBruijnGraph g(k);
            debruijn_graph::FastGraphFromSequencesConstructor<debruijn_graph::DeBruijnGraph>(k, ext_index).ConstructGraph(g, edge_sequences);

            INFO("Saving graph to " << cfg.outfile);
            if (cfg.mode == output_type::gfa) {
                std::ofstream f(cfg.outfile);
                gfa::GFAWriter gfa_writer(g, f);
                gfa_writer.WriteSegmentsAndLinks();
            } else if (cfg.mode == output_type::fastg) {
                io::FastgWriter fastg_writer(g, cfg.outfile);
                fastg_writer.WriteSegmentsAndLinks();
            } else if (cfg.mode == output_type::spades) {
                debruijn_graph::graphio::ConjugateDataPrinter<debruijn_graph::DeBruijnGraph> printer(g);
                debruijn_graph::graphio::PrintBasicGraph(cfg.outfile, printer);
            } else
                FATAL_ERROR("Invalid mode");
        }
    } catch (const std::string &s) {
        std::cerr << s << std::endl;
        return EINTR;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EINTR;
    }

    INFO("SPAdes standalone graph builder finished");

    return 0;
}
