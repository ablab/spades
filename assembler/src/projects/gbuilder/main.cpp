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

#include "io/dataset_support/read_converter.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "io/reads/osequencestream.hpp"

#include "assembly_graph/construction/debruijn_graph_constructor.hpp"

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

struct gcfg {
    gcfg()
            : k(21), outdir("."), nthreads(omp_get_max_threads() / 2 + 1), buff_size(512ULL << 20)
    {}

    unsigned k;
    std::string file;
    std::string outdir;
    unsigned nthreads;
    size_t buff_size;
};


void process_cmdline(int argc, char **argv, gcfg &cfg) {
  using namespace clipp;

  auto cli = (
      cfg.file << value("dataset description (in YAML)"),
      (option("-k") & integer("value", cfg.k)) % "k-mer length to use",
      (option("-t") & integer("value", cfg.nthreads)) % "# of threads to use",
      (option("-o") & value("dir", cfg.outdir)) % "scratch directory to use",
      (option("-b") & integer("value", cfg.buff_size)) % "sorting buffer size, per thread"
  );

  if (!parse(argc, argv, cli)) {
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
        std::string outdir = cfg.outdir, dataset_desc = cfg.file;
        size_t buff_size = cfg.buff_size;

        create_console_logger();

        INFO("Starting SPAdes standalone graph builder, built from " SPADES_GIT_REFSPEC ", git revision " SPADES_GIT_SHA1);

        INFO("K-mer length set to " << k);
        INFO("# of threads to use: " << nthreads);

        nthreads = std::min(nthreads, (unsigned) omp_get_max_threads());
        // Inform OpenMP runtime about this :)
        omp_set_num_threads((int) nthreads);

        fs::make_dir(outdir);
        auto workdir = fs::tmp::make_temp_dir(outdir, "construction");

        io::DataSet<debruijn_graph::config::LibraryData> dataset;
        dataset.load(dataset_desc);
        // FIXME: Get rid of this "/" junk
        debruijn_graph::config::init_libs(dataset, nthreads, buff_size, outdir + "/");

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

        // Step 3: output stuff
        size_t idx = 1;
        for (const auto &edge: edge_sequences) {
            std::cout << std::string(">") << idx++ << "\n";
            io::WriteWrapped(edge.str(), std::cout);
        }
    } catch (const std::string &s) {
        std::cerr << s;
        return EINTR;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << '\n';
        return EINTR;
    }

    return 0;
}
