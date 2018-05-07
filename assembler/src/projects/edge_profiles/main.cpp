//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "io/dataset_support/read_converter.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "io/reads/osequencestream.hpp"

#include "pipeline/graphio.hpp"
#include "pipeline/graph_pack.hpp"

#include "projects/mts/contig_abundance.hpp"
#include "projects/spades/series_analysis.hpp"

#include "utils/stl_utils.hpp"
#include "utils/parallel/openmp_wrapper.h"
#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"

#include <string>
#include <numeric>

#include <cxxopts/cxxopts.hpp>
#include <sys/types.h>
#include <sys/stat.h>

#include "version.hpp"

void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

namespace debruijn_graph {

template<class Graph>
class EdgeProfileStorage : public omnigraph::GraphActionHandler<Graph> {
    typedef omnigraph::GraphActionHandler<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef Profile<Abundance> AbundanceVector;
    
    size_t sample_cnt_;
    std::map<EdgeId, AbundanceVector> profiles_;

    //FIXME self-conjugate edge coverage?!
    template<class SingleStream, class Mapper>
    void Fill(SingleStream &reader, size_t stream_id, const Mapper &mapper) {
        typename SingleStream::ReadT read;
        while (!reader.eof()) {
            reader >> read;
            //TRACE("Aligning " << read.name());

            for (const auto &e_mr: mapper.MapSequence(read.sequence())) {
                // FIXME: should initial_range be used instead?
                profiles_[e_mr.first][stream_id] += double(e_mr.second.mapped_range.size());
            }
        }
    };

public:
    EdgeProfileStorage(const Graph &g, size_t sample_cnt) :
            base(g, "EdgeProfiles"), sample_cnt_(sample_cnt) {}

    template<class SingleStreamList, class Mapper>
    void Fill(SingleStreamList &streams, const Mapper &mapper) {
        for (auto it = this->g().ConstEdgeBegin(); !it.IsEnd(); ++it) {
            profiles_[*it] = AbundanceVector(sample_cnt_, 0.);
        }
        
#       pragma omp parallel for
        for (size_t i = 0; i < sample_cnt_; ++i) {
            Fill(streams[i], i, mapper);
        }
        
        for (auto it = this->g().ConstEdgeBegin(); !it.IsEnd(); ++it) {
            EdgeId e = *it;
            auto &p = profiles_[e];
            for (size_t i = 0; i < sample_cnt_; ++i) {
                p[i] = p[i] / double(this->g().length(e));
            }
        }
    }

    const AbundanceVector& profile(EdgeId e) const {
        return utils::get(profiles_, e);
    }

};

void Run(const string &saves_path, const string &dataset_desc, size_t K,
         const string &profiles_fn, size_t nthreads, const string &tmp_dir) {
    fs::make_dir(tmp_dir);

    io::DataSet<debruijn_graph::config::LibraryData> dataset;
    dataset.load(dataset_desc);

    fs::make_dir(tmp_dir);
    size_t buff_size = 512;
    buff_size <<= 20;
    debruijn_graph::config::init_libs(dataset, nthreads, buff_size, tmp_dir);
    std::vector<size_t> libs(dataset.lib_count());
    std::iota(libs.begin(), libs.end(), 0);

    conj_graph_pack gp(K, tmp_dir, 0);
    gp.kmer_mapper.Attach();
    graphio::ScanGraphPack(saves_path, gp);

    io::BinarySingleStreams single_readers =
            io::single_binary_readers_for_libs(dataset, libs,
                                               /*followed by rc*/true, /*including paired*/true);
    size_t sample_cnt = dataset.lib_count();
    EdgeProfileStorage<Graph> profile_storage(gp.g, sample_cnt);

    profile_storage.Fill(single_readers, *MapperInstance(gp));

    std::ofstream os(profiles_fn);
    for (auto it = gp.g.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        os << gp.g.int_id(e) << PrintVector(profile_storage.profile(e)) << std::endl;
    }
}

}

int main(int argc, char** argv) {
    srand(42);
    srandom(42);
    try {
        unsigned nthreads;
        unsigned k;
        std::string workdir, saves_path, dataset_fn, profiles_fn;

        cxxopts::Options options(argv[0], " analyze strain mix in existing graph");
        options.add_options()
                ("k,kmer", "K-mer length", cxxopts::value<unsigned>(k)->default_value("55"), "K")
                ("g,graph", "saves of the graph", cxxopts::value<std::string>(saves_path))
                ("d,dataset", "Dataset description (in YAML)", cxxopts::value<std::string>(dataset_fn), "file")
                ("p,profiles", "File to store profiles", cxxopts::value<std::string>(profiles_fn), "file")
                ("t,threads", "# of threads to use (default: max threads / 2", cxxopts::value<unsigned>(nthreads)->default_value(std::to_string(omp_get_max_threads() / 2)), "num")
                ("w,workdir", "Working directory (default: .)", cxxopts::value<std::string>(workdir)->default_value("."), "dir")
                ("h,help", "Print help");

        options.parse(argc, argv);
        if (options.count("help")) {
            std::cout << options.help() << std::endl;
            exit(0);
        }

        if (!options.count("graph")) {
            std::cerr << "ERROR: No input files were specified" << std::endl << std::endl;
            std::cout << options.help() << std::endl;
            exit(-1);
        }

        create_console_logger();

        INFO("Built from " SPADES_GIT_REFSPEC ", git revision " SPADES_GIT_SHA1);

        INFO("K-mer length set to " << k);
        INFO("# of threads to use: " << nthreads);

        nthreads = std::min(nthreads, (unsigned) omp_get_max_threads());
        // Inform OpenMP runtime about this :)
        omp_set_num_threads((int) nthreads);

        debruijn_graph::Run(saves_path, dataset_fn, k,
                            profiles_fn, nthreads, workdir + "/tmp/");
    } catch (std::string const &s) {
        std::cerr << s;
        return EINTR;
    } catch (const cxxopts::OptionException &e) {
        std::cerr << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }
}
