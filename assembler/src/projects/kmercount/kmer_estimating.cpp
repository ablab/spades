//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "adt/cyclichash.hpp"

#include "version.hpp"
#include "io/dataset_support/read_converter.hpp"
#include "io/reads/osequencestream.hpp"
#include "io/reads/coverage_filtering_read_wrapper.hpp"
#include "utils/parallel/openmp_wrapper.h"
#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"
#include "utils/kmer_counting.hpp"

#include <sys/types.h>
#include <string>
#include <clipp/clipp.h>

using namespace std;
void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

namespace kmer_estimating {
struct Args {
    unsigned k = 21;
    std::string dataset_desc;
    unsigned nthreads = (omp_get_max_threads() / 2);
};
}

void process_cmdline(int argc, char **argv, kmer_estimating::Args &args) {
    using namespace clipp;
    bool print_help = false;

    auto cli = (
        (option("-k", "--kmer") & integer("value", args.k)) % "K-mer length",
        (required("-d", "--dataset") & value("dir", args.dataset_desc)) % "Dataset description (in YAML)",
        (option("-t", "--threads") & integer("value", args.nthreads)) % "# of threads to use",
        (option("-h", "--help").set(print_help)) % "Show help"
    );

    auto result = parse(argc, argv, cli);
    if (!result || print_help) {
        std::cout << make_man_page(cli, argv[0]).prepend_section("DESCRIPTION", " Kmer number estimating.  Kmers from reverse-complementary reads aren't taken into account.");
        if (print_help) {
            exit(0);
        } else {
            exit(1);
        }
    }
}

int main(int argc, char* argv[]) {
    utils::perf_counter pc;

    srand(42);
    srandom(42);
    try {
        std::vector<std::string> input;

        kmer_estimating::Args args;
        process_cmdline(argc, argv, args);

        create_console_logger();

        args.nthreads = spades_set_omp_threads(args.nthreads);

        INFO("Starting kmer spectra cardinality, built from " << version::refspec() << ", git revision " << version::gitrev());

        INFO("K-mer length set to " << args.k);
        INFO("# of threads to use: " << args.nthreads);
        INFO("Maximum # of threads to use (adjusted due to OMP capabilities): " << args.nthreads);

        io::DataSet<debruijn_graph::config::LibraryData> dataset;
        dataset.load(args.dataset_desc);

        std::vector<size_t> libs(dataset.lib_count());
        std::iota(libs.begin(), libs.end(), 0);

        io::SingleStreams single_readers;

        for (size_t i = 0; i < dataset.lib_count(); ++i) {
            single_readers.push_back(io::single_easy_reader(dataset[i], false, true));
        }

        INFO("Estimating kmer cardinality");
        typedef rolling_hash::SymmetricCyclicHash<> SeqHasher;
        SeqHasher hasher(args.k);

        size_t kmers_cnt_est = utils::EstimateCardinalityForOneStream(args.k, single_readers, hasher);
        INFO("Kmer number estimation: " << kmers_cnt_est)
    } catch (std::string const &s) {
        std::cerr << s;
        return EINTR;
    }

    return 0;
}
