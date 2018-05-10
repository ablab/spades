//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <cxxopts/cxxopts.hpp>
//#include <omp.h>
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

template<class IS, class OS>
void Transfer(IS &input, OS &output) {
    typename OS::ReadT read;
    while (!input.eof()) {
        input >> read;
        output << read;
    }
}

int main(int argc, char* argv[]) {
    typedef qf::cqf CQFKmerFilter;
    //typedef CyclicHash<64, uint8_t, NDNASeqHash<uint8_t>> SeqHasher;
    utils::perf_counter pc;

    srand(42);
    srandom(42);
    try {
        unsigned nthreads;
        unsigned thr, k;
        std::string workdir, dataset_desc;
        std::vector<std::string> input;
        size_t buff_size = 512;
        buff_size <<= 20;

        cxxopts::Options options(argv[0], " kmer count read filter");
        options.add_options()
                ("k,kmer", "K-mer length", cxxopts::value<unsigned>(k)->default_value("21"), "K")
                ("c,cov", "Median kmer count threshold (read pairs, s.t. kmer count median for BOTH reads LESS OR EQUAL to this value will be ignored)", cxxopts::value<unsigned>(thr)->default_value("2"), "threshold")
                ("d,dataset", "Dataset description (in YAML)", cxxopts::value<std::string>(dataset_desc), "file")
                ("t,threads", "# of threads to use", cxxopts::value<unsigned>(nthreads)->default_value(std::to_string(omp_get_max_threads() / 2)), "num")
                ("o,outdir", "Output directory to use", cxxopts::value<std::string>(workdir)->default_value("."), "dir")
//                ("b,bufsize", "Sorting buffer size, per thread", cxxopts::value<size_t>(read_buffer_size)->default_value("536870912"))
                ("h,help", "Print help");

//        options.add_options("Input")
//                ("positional", "", cxxopts::value<std::vector<std::string>>(input));
//
//        options.parse_positional("positional");
        options.parse(argc, argv);
        if (options.count("help")) {
            std::cout << options.help() << std::endl;
            exit(0);
        }

//        if (!options.count("positional") && !options.count("dataset")) {
//            std::cerr << "ERROR: No input files were specified" << std::endl << std::endl;
//            std::cout << options.help() << std::endl;
//            exit(-1);
//        }
        if (!options.count("dataset")) {
            std::cerr << "ERROR: No input files were specified" << std::endl << std::endl;
            std::cout << options.help() << std::endl;
            exit(-1);
        }

        create_console_logger();

        nthreads = spades_set_omp_threads(nthreads);

        INFO("Starting kmer count based read filtering, built from " SPADES_GIT_REFSPEC ", git revision " SPADES_GIT_SHA1);

        INFO("K-mer length set to " << k);
        INFO("# of threads to use: " << nthreads);
        INFO("Maximum # of threads to use (adjusted due to OMP capabilities): " << nthreads);

        io::DataSet<debruijn_graph::config::LibraryData> dataset;
        dataset.load(dataset_desc);

        fs::make_dirs(workdir + "/tmp/");
        debruijn_graph::config::init_libs(dataset, nthreads, buff_size, workdir + "/tmp/");

        std::vector<size_t> libs(dataset.lib_count());
        std::iota(libs.begin(), libs.end(), 0);
        io::BinarySingleStreams single_readers = io::single_binary_readers_for_libs(dataset, libs,
                                                                                    /*followed by rc*/false, /*including paired*/true);
        INFO("Estimating kmer cardinality");
        typedef rolling_hash::SymmetricCyclicHash<> SeqHasher;
        SeqHasher hasher(k);

        size_t kmers_cnt_est = utils::EstimateCardinality(k, single_readers, hasher);

        CQFKmerFilter cqf(kmers_cnt_est);
        INFO("Filling kmer coverage");
        utils::FillCoverageHistogram(cqf, k, hasher, single_readers, thr + 1);
        INFO("Kmer coverage filled");

        for (size_t i = 0; i < dataset.lib_count(); ++i) {
            INFO("Filtering library " << i);
            if (dataset[i].has_paired()) {
                auto filtered = io::CovFilteringWrap<io::PairedRead>(
                        io::paired_easy_reader(dataset[i], /*followed by rc*/false, /*insert size*/0),
                        k, hasher, cqf, thr);
                io::OFastqPairedStream ostream(workdir + "/" + to_string(i + 1) + ".1.fastq",
                                              workdir + "/" + to_string(i + 1) + ".2.fastq");
                Transfer(*filtered, ostream);
            }

            if (dataset[i].has_single()) {
                auto filtered = io::CovFilteringWrap<io::SingleRead>(
                        io::single_easy_reader(dataset[i], /*followed_by_rc*/ false,
                                               /*including_paired_reads*/ false),
                        k, hasher, cqf, thr);

                io::OFastqReadStream ostream(workdir + "/" + to_string(i + 1) + ".s.fastq");
                Transfer(*filtered, ostream);
            }
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
