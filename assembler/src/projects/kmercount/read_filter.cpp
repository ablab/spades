//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <sys/types.h>
#include <string>
#include <cxxopts/cxxopts.hpp>

#include "utils/parallel/openmp_wrapper.h"
#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"
#include "io/dataset_support/read_converter.hpp"
#include "adt/cyclichash.hpp"
#include "adt/cqf.hpp"
#include "io/reads/osequencestream.hpp"
#include "utils/kmer_counting.hpp"
#include "io/reads/coverage_filtering_read_wrapper.hpp"

#include "version.hpp"

using namespace std;
typedef rolling_hash::SymmetricCyclicHash<> SeqHasher;

void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

template<class IS, class OS, class Filter>
void filter_reads(IS &input, OS &output, const Filter& filter, unsigned buffer_size, unsigned nthreads) {
    std::vector<typename OS::ReadT> reads_buffer(buffer_size);
    std::vector<uint8_t> need_to_out(buffer_size);
    std::vector<unsigned> chunk_start(nthreads), chunk_end(nthreads);

    while (!input.eof()) {
        unsigned reads_cnt = 0;
        while (!input.eof() && reads_cnt < reads_buffer.size()) {
            input >> reads_buffer[reads_cnt];
            ++reads_cnt;
        }

        unsigned reads_per_thread = reads_cnt/nthreads;
        chunk_start[0] = 0;
        chunk_end[0] = reads_per_thread;
        for (unsigned i = 1; i < nthreads; ++i) {
            chunk_start[i] = chunk_end[i - 1];
            chunk_end[i] = chunk_start[i] + reads_per_thread;
        }
        chunk_end[nthreads - 1] = reads_cnt;

#       pragma omp parallel for
        for (unsigned i = 0; i < nthreads; ++i) {
            for (unsigned j = chunk_start[i]; j < chunk_end[i]; ++j) {
                if (!filter(reads_buffer[j])) {
                    need_to_out[j] = 1;
                } else {
                    need_to_out[j] = 0;
                }
            }
        }

        for (size_t i = 0; i < reads_cnt; ++i) {
            if (need_to_out[i] == 1) {
                output << reads_buffer[i];
            }
            need_to_out[i] = 0;
        }
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

        INFO("Starting kmer count based read filtering, built from " << version::refspec() << ", git revision " << version::gitrev());

        INFO("K-mer length set to " << k);
        INFO("# of threads to use: " << nthreads);
        INFO("Maximum # of threads to use (adjusted due to OMP capabilities): " << nthreads);

        io::DataSet<debruijn_graph::config::LibraryData> dataset;
        dataset.load(dataset_desc);

        fs::make_dirs(workdir + "/tmp/");
        debruijn_graph::config::init_libs(dataset, nthreads, buff_size, workdir + "/tmp/");

        for (size_t i = 0; i < dataset.lib_count(); ++i) {
            io::ReadConverter::ConvertToBinary(dataset[i]);
        }

        std::vector<size_t> libs(dataset.lib_count());
        std::iota(libs.begin(), libs.end(), 0);
        io::BinarySingleStreams single_readers = io::single_binary_readers_for_libs(dataset, libs,
                                                                                    /*followed by rc*/false, /*including paired*/true);
        INFO("Estimating kmer cardinality");
        SeqHasher hasher(k);

        size_t kmers_cnt_est = utils::EstimateCardinality(k, single_readers, hasher);

        CQFKmerFilter cqf(kmers_cnt_est);
        INFO("Filling kmer coverage");
        utils::FillCoverageHistogram(cqf, k, hasher, single_readers, thr + 1);
        INFO("Kmer coverage filled");

        const unsigned FILTER_READS_BUFF_SIZE = 1 << 20;

        for (size_t i = 0; i < dataset.lib_count(); ++i) {
            INFO("Filtering library " << i);
            if (dataset[i].has_paired()) {
                io::PairedStreamPtr paired_reads_stream =
                    io::paired_easy_reader(dataset[i], /*followed by rc*/false, /*insert size*/0);
                io::OFastqPairedStream ostream(workdir + "/" + to_string(i + 1) + ".1.fastq",
                                               workdir + "/" + to_string(i + 1) + ".2.fastq");
                io::CoverageFilter<io::PairedRead, SeqHasher> filter(k, hasher, cqf, thr);
                filter_reads(*paired_reads_stream, ostream, filter, FILTER_READS_BUFF_SIZE, nthreads);
            }

            if (dataset[i].has_single()) {
                io::SingleStreamPtr single_reads_stream = io::single_easy_reader(dataset[i],
                    /*followed_by_rc*/ false, /*including_paired_reads*/ false);
                io::CoverageFilter<io::SingleRead, SeqHasher> filter(k, hasher, cqf, thr);

                io::OFastqReadStream ostream(workdir + "/" + to_string(i + 1) + ".s.fastq");
                filter_reads(*single_reads_stream, ostream, filter, FILTER_READS_BUFF_SIZE, nthreads);
            }
        }
        INFO("Filtering finished")
    } catch (std::string const &s) {
        std::cerr << s;
        return EINTR;
    } catch (const cxxopts::OptionException &e) {
        std::cerr << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }

    return 0;
}
