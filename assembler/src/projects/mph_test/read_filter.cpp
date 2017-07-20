//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <cxxopts/cxxopts.hpp>
#include <omp.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string>
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
    typedef qf::cqf<RtSeq> CQFKmerFilter;
    typedef CyclicHash<64, uint8_t, NDNASeqHash<uint8_t>> SeqHasher;
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
                ("c,cov", "Median kmer count threshold", cxxopts::value<unsigned>(thr)->default_value("2"), "threshold")
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

        INFO("Starting kmer count based read filtering, built from " SPADES_GIT_REFSPEC ", git revision " SPADES_GIT_SHA1);

        INFO("K-mer length set to " << k);
        INFO("# of threads to use: " << nthreads);

        nthreads = std::min(nthreads, (unsigned) omp_get_max_threads());
        // Inform OpenMP runtime about this :)
        omp_set_num_threads((int) nthreads);

        io::DataSet<debruijn_graph::config::DataSetData> dataset;
        dataset.load(dataset_desc);

        fs::make_dir(workdir + "/tmp/");
        debruijn_graph::config::init_libs(dataset, nthreads, buff_size, workdir + "/tmp/");

        utils::StoringTypeFilter<utils::InvertableStoring> filter;
        SeqHasher hasher(k);
        io::BinarySingleStreams single_readers = io::single_binary_readers(dataset, /*followed by rc*/false, /*including paired*/true);
        size_t kmers_cnt_est = EstimateCardinality(k, single_readers, filter);
        CQFKmerFilter cqf([=](const RtSeq &s) { return hasher.hash(s); },
                                              kmers_cnt_est);

        utils::FillCoverageHistogram(cqf, k, single_readers, filter, thr + 1);
        
        auto filter_f = [=,&cqf] (io::PairedRead& p_r) { return io::CountMedianMlt(p_r.first().sequence(), k, hasher, cqf) > thr ||
                                                           io::CountMedianMlt(p_r.second().sequence(), k, hasher, cqf) > thr; };

        for (size_t i = 0; i < dataset.lib_count(); ++i) {
            auto filtered = io::FilteringWrap<io::PairedRead>(io::paired_easy_library_reader(dataset[i], 
                        /*followed by rc*/false, /*insert size*/0),
                                                              filter_f);
            io::OPairedReadStream ostream(workdir + "/" + to_string(i + 1) + ".1.fastq",
                                          workdir + "/" + to_string(i + 1) + ".2.fastq");

            //FIXME find utility method?
            io::PairedRead paired_read;
            while (!filtered->eof()) {
                (*filtered) >> paired_read;
                ostream << paired_read;
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
