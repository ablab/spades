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

#include "assembly_graph/construction/debruijn_graph_constructor.hpp"

#include "version.hpp"

#include <cxxopts/cxxopts.hpp>

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

int main(int argc, char* argv[]) {
    srand(42);
    srandom(42);

    try {
        unsigned nthreads;
        unsigned k;
        std::string outdir, dataset_desc;
        std::vector<std::string> input;
        size_t buff_size = 512;
        buff_size <<= 20;

        cxxopts::Options options(argv[0], " kmer count read filter");
        options.add_options()
                ("k,kmer", "K-mer length", cxxopts::value<unsigned>(k)->default_value("21"), "K")
                ("d,dataset", "Dataset description (in YAML)", cxxopts::value<std::string>(dataset_desc), "file")
                ("t,threads", "# of threads to use", cxxopts::value<unsigned>(nthreads)->default_value(std::to_string(omp_get_max_threads() / 2 + 1)), "num")
                ("o,outdir", "Output directory to use", cxxopts::value<std::string>(outdir)->default_value("."), "dir")
                ("b,bufsize", "Sorting buffer size, per thread", cxxopts::value<size_t>(buff_size)->default_value("536870912"))
                ("h,help", "Print help");

        options.parse(argc, argv);
        if (options.count("help")) {
            std::cout << options.help() << std::endl;
            exit(0);
        }

        if (!options.count("dataset")) {
            std::cerr << "ERROR: No input files were specified" << std::endl << std::endl;
            std::cout << options.help() << std::endl;
            exit(-1);
        }

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
            std::cerr << ">" << idx++ << '\n';
            std::cerr << edge << '\n';
        }
    } catch (std::string const &s) {
        std::cerr << s;
        return EINTR;
    } catch (const cxxopts::OptionException &e) {
        std::cerr << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }

    return 0;
}
