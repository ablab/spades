//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "version.hpp"

#include "adt/cyclichash.hpp"
#include "adt/cqf.hpp"

#include "io/dataset_support/read_converter.hpp"
#include "io/reads/osequencestream.hpp"
#include "io/reads/coverage_filtering_read_wrapper.hpp"

#include "utils/parallel/openmp_wrapper.h"
#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"
#include "utils/kmer_counting.hpp"

#include "threadpool/threadpool.hpp"

#include <clipp/clipp.h>
#include <sys/types.h>
#include <string>

using namespace std;
typedef rolling_hash::SymmetricCyclicHash<> SeqHasher;

void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

namespace read_filter { 
struct Args {
    unsigned thr = 2, k = 21;
    std::string dataset_desc, workdir = ".";
    unsigned nthreads = (omp_get_max_threads() / 2);
};
}

void process_cmdline(int argc, char **argv, read_filter::Args &args) {
    using namespace clipp;
    bool print_help = false;
    
    auto cli = (
        (option("-k", "--kmer") & integer("value", args.k)) % "K-mer length",
        (option("-c", "--cov") & integer("value", args.thr)) % "Median kmer count threshold (read pairs, s.t. kmer count median for BOTH reads LESS OR EQUAL to this value will be ignored)",
        (required("-d", "--dataset") & value("yaml", args.dataset_desc)) % "Dataset description (in YAML)",
        (option("-t", "--threads") & integer("value", args.nthreads)) % "# of threads to use",
        (option("-o", "--outdir") & value("dir", args.workdir)) %  "Output directory to use",
        (option("-h", "--help").set(print_help)) % "Show help"
    );

    auto result = parse(argc, argv, cli);
    if (!result || print_help) {
        std::cout << make_man_page(cli, argv[0]).prepend_section("DESCRIPTION", " Kmer count read filter");
        if (print_help) {
            exit(0);
        } else {
            exit(1);
        }
    }
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
                if (filter(reads_buffer[j])) {
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
        std::vector<std::string> input;

        read_filter::Args args;
        process_cmdline(argc, argv, args);

        create_console_logger();

        args.nthreads = spades_set_omp_threads(args.nthreads);

        INFO("Starting kmer count based read filtering, built from " << version::refspec() << ", git revision " << version::gitrev());

        INFO("K-mer length set to " << args.k);
        INFO("# of threads to use: " << args.nthreads);
        INFO("Maximum # of threads to use (adjusted due to OMP capabilities): " << args.nthreads);

        io::DataSet<debruijn_graph::config::LibraryData> dataset;
        dataset.load(args.dataset_desc);

        fs::make_dirs(args.workdir + "/tmp/");
        debruijn_graph::config::init_libs(dataset, args.nthreads, args.workdir + "/tmp/");

        std::unique_ptr<ThreadPool::ThreadPool> pool;

        if (args.nthreads > 1) {
            pool = std::make_unique<ThreadPool::ThreadPool>(args.nthreads);
        }

        for (size_t i = 0; i < dataset.lib_count(); ++i) {
            io::ReadConverter::ConvertToBinary(dataset[i], pool.get());
        }

        std::vector<size_t> libs(dataset.lib_count());
        std::iota(libs.begin(), libs.end(), 0);
        io::BinarySingleStreams single_readers = io::single_binary_readers_for_libs(dataset, libs,
                                                                                    /*followed by rc*/false, /*including paired*/true);
        INFO("Estimating kmer cardinality");
        SeqHasher hasher(args.k);

        size_t kmers_cnt_est = utils::EstimateCardinalityUpperBound(args.k, single_readers, hasher);

        CQFKmerFilter cqf(kmers_cnt_est);
        INFO("Filling kmer coverage");
        utils::FillCoverageHistogram(cqf, args.k, hasher, single_readers, args.thr + 1);
        INFO("Kmer coverage filled");

        const unsigned FILTER_READS_BUFF_SIZE = 1 << 20;

        for (size_t i = 0; i < dataset.lib_count(); ++i) {
            INFO("Filtering library " << i);
            if (dataset[i].has_paired()) {
                io::PairedStream paired_reads_stream =
                        io::paired_easy_reader(dataset[i], /*followed by rc*/false, /*insert size*/0);
                io::OFastqPairedStream ostream(args.workdir + "/" + to_string(i + 1) + ".1.fastq",
                                               args.workdir + "/" + to_string(i + 1) + ".2.fastq");
                io::CoverageFilter<io::PairedRead, SeqHasher> filter(args.k, hasher, cqf, args.thr + 1);
                filter_reads(paired_reads_stream, ostream, filter, FILTER_READS_BUFF_SIZE, args.nthreads);
            }

            if (dataset[i].has_single()) {
                io::SingleStream single_reads_stream =
                        io::single_easy_reader(dataset[i], /*followed_by_rc*/ false, /*including_paired_reads*/ false);
                io::CoverageFilter<io::SingleRead, SeqHasher> filter(args.k, hasher, cqf, args.thr + 1);

                io::OFastqReadStream ostream(args.workdir + "/" + to_string(i + 1) + ".s.fastq");
                filter_reads(single_reads_stream, ostream, filter, FILTER_READS_BUFF_SIZE, args.nthreads);
            }
        }
        INFO("Filtering finished")
    } catch (std::string const &s) {
        std::cerr << s;
        return EINTR;
    }

    return 0;
}
