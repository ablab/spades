//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <cxxopts/cxxopts.hpp>
#include <sys/types.h>
#include <sys/stat.h>
#include <string>
#include "utils/parallel/openmp_wrapper.h"
#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"

#include "io/dataset_support/read_converter.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "pipeline/config_struct.hpp"
#include "adt/cyclichash.hpp"
#include "adt/cqf.hpp"
#include "io/reads/osequencestream.hpp"
#include "utils/kmer_counting.hpp"
#include "utils/ph_map/storing_traits.hpp"
#include "io/reads/coverage_filtering_read_wrapper.hpp"

#include "version.hpp"

using namespace std;
void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main(int argc, char* argv[]) {
    utils::perf_counter pc;

    srand(42);
    srandom(42);
    try {
        unsigned nthreads;
        unsigned k;
        std::string dataset_desc;
        std::vector<std::string> input;
        size_t buff_size = 512;
        buff_size <<= 20;

        cxxopts::Options options(argv[0], " kmer number estimating.  Kmers from reverse-complementary reads aren't taken into account.");
        options.add_options()
            ("k,kmer", "K-mer length", cxxopts::value<unsigned>(k)->default_value("21"), "K")
            ("d,dataset", "Dataset description (in YAML)", cxxopts::value<std::string>(dataset_desc), "file")
            ("t,threads", "# of threads to use", cxxopts::value<unsigned>(nthreads)->default_value(std::to_string(omp_get_max_threads() / 2)), "num")
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

        nthreads = spades_set_omp_threads(nthreads);

        INFO("Starting kmer spectra cardinality, built from " << version::refspec() << ", git revision " << version::gitrev());

        INFO("K-mer length set to " << k);
        INFO("# of threads to use: " << nthreads);
        INFO("Maximum # of threads to use (adjusted due to OMP capabilities): " << nthreads);

        io::DataSet<debruijn_graph::config::LibraryData> dataset;
        dataset.load(dataset_desc);

        std::vector<size_t> libs(dataset.lib_count());
        std::iota(libs.begin(), libs.end(), 0);

        io::SingleStreams single_readers;

        for (size_t i = 0; i < dataset.lib_count(); ++i) {
            single_readers.push_back(io::single_easy_reader(dataset[i], false, true));
        }

        INFO("Estimating kmer cardinality");
        typedef rolling_hash::SymmetricCyclicHash<> SeqHasher;
        SeqHasher hasher(k);

        size_t kmers_cnt_est = utils::EstimateCardinality(k, single_readers, hasher);
        INFO("Kmer number estimation: " << kmers_cnt_est)
    } catch (std::string const &s) {
        std::cerr << s;
        return EINTR;
    } catch (const cxxopts::OptionException &e) {
        std::cerr << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }

    return 0;
}
