//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "io/dataset_support/read_converter.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "io/graph/gfa_reader.hpp"

#include "pipeline/graphio.hpp"
#include "pipeline/graph_pack.hpp"

#include "projects/mts/contig_abundance.hpp"
#include "projects/spades/series_analysis.hpp"

#include "utils/parallel/openmp_wrapper.h"
#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"

#include <clipp/clipp.h>

#include <string>
#include <numeric>

#include <sys/types.h>
#include <sys/stat.h>

#include "version.hpp"

void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

static bool ends_with(const std::string &s, const std::string &p) {
    if (s.size() < p.size())
        return false;

    return (s.compare(s.size() - p.size(), p.size(), p) == 0);
}

void LoadGraph(debruijn_graph::ConjugateDeBruijnGraph &graph, const std::string &filename) {
    using namespace debruijn_graph;
    if (ends_with(filename, ".gfa")) {
        gfa::GFAReader gfa(filename);
        INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links());
        gfa.to_graph(graph);
    } else {
        graphio::ScanBasicGraph(filename, graph);
    }
}

using namespace debruijn_graph;

class EdgeProfileStorage {
    typedef ConjugateDeBruijnGraph::EdgeId EdgeId;
    typedef Profile<Abundance> AbundanceVector;

    const ConjugateDeBruijnGraph &g_;
    size_t sample_cnt_;
    std::map<EdgeId, AbundanceVector> profiles_;

    // FIXME self-conjugate edge coverage?!
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
    EdgeProfileStorage(const ConjugateDeBruijnGraph &g, size_t sample_cnt) :
            g_(g), sample_cnt_(sample_cnt) {}

    template<class SingleStreamList, class Mapper>
    void Fill(SingleStreamList &streams, const Mapper &mapper) {
        for (auto it = g_.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            profiles_[*it] = AbundanceVector(sample_cnt_, 0.);
        }

#       pragma omp parallel for
        for (size_t i = 0; i < sample_cnt_; ++i) {
            Fill(streams[i], i, mapper);
        }

        for (auto it = g_.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            EdgeId e = *it;
            auto &p = profiles_[e];
            for (size_t i = 0; i < sample_cnt_; ++i) {
                p[i] = p[i] / float(g_.length(e));
            }
        }
    }

    const AbundanceVector& profile(EdgeId e) const {
        return utils::get(profiles_, e);
    }

};

typedef io::DataSet<debruijn_graph::config::LibraryData> DataSet;
typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> SequencingLib;

void Run(const std::string &graph_path, const std::string &dataset_desc, size_t K,
         const std::string &profiles_fn, size_t nthreads, const std::string &tmpdir) {
    DataSet dataset;
    dataset.load(dataset_desc);

    debruijn_graph::conj_graph_pack gp(K, tmpdir, dataset.lib_count());

    INFO("Loading de Bruijn graph from " << graph_path);
    LoadGraph(gp.g, graph_path);
    {
        size_t sz = 0;
        for (auto it = gp.g.ConstEdgeBegin(); !it.IsEnd(); ++it)
            sz += 1;

        INFO("Graph loaded. Total vertices: " << gp.g.size() << " Total edges: " << sz);
    }

    // FIXME: Get rid of this "/" junk
    debruijn_graph::config::init_libs(dataset, nthreads, 512ULL << 20, tmpdir + "/");

    gp.kmer_mapper.Attach();
    gp.EnsureBasicMapping();

    std::vector<size_t> libs(dataset.lib_count());
    std::iota(libs.begin(), libs.end(), 0);

    io::BinarySingleStreams single_readers =
            io::single_binary_readers_for_libs(dataset, libs,
                                               /*followed by rc*/true, /*including paired*/true);
    size_t sample_cnt = dataset.lib_count();
    EdgeProfileStorage profile_storage(gp.g, sample_cnt);

    profile_storage.Fill(single_readers, *MapperInstance(gp));

    std::ofstream os(profiles_fn);
    for (auto it = gp.g.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        os << gp.g.int_id(e) << '\t' << PrintVector(profile_storage.profile(e)) << std::endl;
    }
}

struct gcfg {
    gcfg()
        : k(21), tmpdir("tmp"), outfile("-"),
          nthreads(omp_get_max_threads() / 2 + 1)
    {}

    unsigned k;
    std::string file;
    std::string graph;
    std::string tmpdir;
    std::string outfile;
    unsigned nthreads;
};

void process_cmdline(int argc, char **argv, gcfg &cfg) {
  using namespace clipp;

  auto cli = (
      cfg.file << value("dataset description (in YAML)"),
      cfg.graph << value("graph (in GFA)"),
      cfg.outfile << value("output filename"),
      (option("-k") & integer("value", cfg.k)) % "k-mer length to use",
      (option("-t") & integer("value", cfg.nthreads)) % "# of threads to use",
      (option("-tmpdir") & value("dir", cfg.tmpdir)) % "scratch directory to use"
  );

  auto result = parse(argc, argv, cli);
  if (!result) {
      std::cout << make_man_page(cli, argv[0]);
      exit(1);
  }
}

int main(int argc, char** argv) {
    utils::segfault_handler sh;
    gcfg cfg;

    process_cmdline(argc, argv, cfg);

    create_console_logger();
    INFO("Built from " SPADES_GIT_REFSPEC ", git revision " SPADES_GIT_SHA1);

    try {
        unsigned nthreads = cfg.nthreads;
        unsigned k = cfg.k;
        std::string tmpdir = cfg.tmpdir;

        fs::make_dir(tmpdir);

        INFO("K-mer length set to " << k);
        INFO("# of threads to use: " << nthreads);

        nthreads = std::min(nthreads, (unsigned) omp_get_max_threads());
        // Inform OpenMP runtime about this :)
        omp_set_num_threads((int) nthreads);

        Run(cfg.graph, cfg.file, k, cfg.outfile, nthreads, tmpdir);
    } catch (const std::string &s) {
        std::cerr << s << std::endl;
        return EINTR;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EINTR;
    }
}
