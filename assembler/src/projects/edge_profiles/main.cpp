//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "profile_storage.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "io/graph/gfa_reader.hpp"

#include "io/binary/graph_pack.hpp"
#include "projects/mts/contig_abundance.hpp"

#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"

#include "version.hpp"

#include <clipp/clipp.h>
#include <unordered_map>
#include <string>
#include <numeric>
#include <sys/types.h>
#include <sys/stat.h>

static void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

using namespace debruijn_graph;

static void PrintGraphInfo(Graph &g) {
    size_t sz = 0;
    for (auto it = g.ConstEdgeBegin(); !it.IsEnd(); ++it)
        sz += 1;

    INFO("Graph loaded. Total vertices: " << g.size() << " Total edges: " << sz);
}

static void LoadGraph(conj_graph_pack &gp, const std::string &filename) {
    //if (ends_with(filename, ".gfa")) {
    //    gfa::GFAReader gfa(filename);
    //    INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links());
    //    gfa.to_graph(gp.g);
    //} else {
    io::binary::BasePackIO<Graph>().Load(filename, gp);
    //}
    PrintGraphInfo(gp.g);
}


typedef io::DataSet<config::LibraryData> DataSet;
typedef io::SequencingLibrary<config::LibraryData> SequencingLib;

static io::ReadStreamList<io::SingleRead>
single_easy_readers_for_libs(DataSet& dataset_info,
                             const std::vector<size_t>& libs,
                             bool followed_by_rc = true,
                             bool including_paired_reads = true,
                             bool handle_Ns = true,
                             io::OffsetType offset_type = io::PhredOffset) {
    VERIFY(!libs.empty());
    io::ReadStreamList<io::SingleRead> streams;
    for (auto l_id : libs) {
        streams.push_back(io::single_easy_reader(dataset_info[l_id],
                                             followed_by_rc,
                                             including_paired_reads, handle_Ns, offset_type));
    }
    return streams;
}

static void Run(const std::string &graph_path, const std::string &dataset_desc, size_t K,
         const std::string &profiles_fn, size_t nthreads, const std::string &tmpdir) {
    DataSet dataset;
    dataset.load(dataset_desc);

    conj_graph_pack gp(K, tmpdir, dataset.lib_count());

    INFO("Loading de Bruijn graph from " << graph_path);
    gp.kmer_mapper.Attach();

    LoadGraph(gp, graph_path);

    // FIXME: Get rid of this "/" junk
    config::init_libs(dataset, nthreads, tmpdir + "/");

    gp.EnsureBasicMapping();

    std::vector<size_t> libs(dataset.lib_count());
    std::iota(libs.begin(), libs.end(), 0);

    auto single_readers = single_easy_readers_for_libs(dataset, libs,
            /*followed by rc*/true, /*including paired*/true);

    size_t sample_cnt = dataset.lib_count();
    debruijn_graph::coverage_profiles::EdgeProfileStorage profile_storage(gp.g, sample_cnt);

    profile_storage.Fill(single_readers, *MapperInstance(gp));

    std::ofstream os(profiles_fn);
    profile_storage.Save(os);
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

static void process_cmdline(int argc, char **argv, gcfg &cfg) {
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
    START_BANNER("edge profile counter");

    try {
        unsigned nthreads = cfg.nthreads;
        unsigned k = cfg.k;
        std::string tmpdir = cfg.tmpdir;

        fs::make_dir(tmpdir);

        INFO("K-mer length set to " << k);

        nthreads = std::min(nthreads, (unsigned) omp_get_max_threads());
        // Inform OpenMP runtime about this :)
        omp_set_num_threads((int) nthreads);
        INFO("# of threads to use: " << nthreads);

        Run(cfg.graph, cfg.file, k, cfg.outfile, nthreads, tmpdir);
    } catch (const std::string &s) {
        std::cerr << s << std::endl;
        return EINTR;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EINTR;
    }
}
