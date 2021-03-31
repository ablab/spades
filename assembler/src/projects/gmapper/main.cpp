//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"

#include "modules/alignment/bwa_index.hpp"
#include "modules/alignment/sequence_mapper_notifier.hpp"
#include "modules/alignment/bwa_sequence_mapper.hpp"
#include "modules/alignment/long_read_mapper.hpp"

#include "io/graph/gfa_reader.hpp"
#include "io/dataset_support/read_converter.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "io/binary/graph.hpp"

#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"

#include "paired_info/paired_info_utils.hpp"
#include "utils/logger/log_writers.hpp"
#include "utils/logger/logger.hpp"
#include "utils/segfault_handler.hpp"
#include "utils/filesystem/temporary.hpp"

#include "pipeline/graph_pack.hpp"
#include "pipeline/configs/aligner_config.hpp"

#include "projects/spades/hybrid_aligning.hpp"
#include "projects/spades/hybrid_gap_closer.hpp"

#include "version.hpp"

#include "threadpool/threadpool.hpp"
#include <clipp/clipp.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <limits>
#include <string>
#include <fstream>

void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

using namespace debruijn_graph;

struct gcfg {
    gcfg()
        : k(21), tmpdir("tmp"), outfile("-"),
          nthreads(omp_get_max_threads() / 2 + 1),
          libindex(0),
          mode(alignment::BWAIndex::AlignmentMode::Default)
    {}

    unsigned k;
    std::string file;
    std::string graph;
    std::string tmpdir;
    std::string outfile;
    unsigned nthreads;
    unsigned libindex;
    alignment::BWAIndex::AlignmentMode mode;
};

void process_cmdline(int argc, char **argv, gcfg &cfg) {
  using namespace clipp;

  auto cli = (
      cfg.file << value("dataset description (in YAML)"),
      cfg.graph << value("graph (in GFA)"),
      cfg.outfile << value("output filename"),
      (option("-l") & integer("value", cfg.libindex)) % "library index (0-based, default: 0)",
      (option("-k") & integer("value", cfg.k)) % "k-mer length to use",
      (option("-t") & integer("value", cfg.nthreads)) % "# of threads to use",
      (option("--tmp-dir") & value("dir", cfg.tmpdir)) % "scratch directory to use",
      (with_prefix("-X", option("pacbio").set(cfg.mode, alignment::BWAIndex::AlignmentMode::PacBio) |
                         option("intractg").set(cfg.mode, alignment::BWAIndex::AlignmentMode::IntraCtg)) % "inner alignment mode (for paired-end / contigs)")
  );

  auto result = parse(argc, argv, cli);
  if (!result) {
      auto fmt = doc_formatting{}
                 .merge_alternative_flags_with_common_prefix(true);

      std::cout << make_man_page(cli, argv[0], fmt);
      exit(1);
  }
}

typedef io::DataSet<debruijn_graph::config::LibraryData> DataSet;
typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> SequencingLib;

static void ProcessContigs(const Graph &graph, alignment::BWAIndex::AlignmentMode mode,
                           SequencingLib &lib,
                           PathStorage<Graph> &single_long_reads,
                           path_extend::GappedPathStorage &trusted_paths) {
    SequenceMapperNotifier notifier;
    LongReadMapper read_mapper(graph, single_long_reads, trusted_paths, lib.type());
    notifier.Subscribe(&read_mapper);

    INFO("Mapping using BWA-mem mapper");
    alignment::BWAReadMapper<Graph> mapper(graph, mode);
    auto single_streams = single_easy_readers(lib, false, false, /*handle Ns*/ false);
    notifier.ProcessLibrary(single_streams, mapper);
}

void LoadGraph(debruijn_graph::ConjugateDeBruijnGraph &graph, const std::string &filename,
               io::IdMapper<std::string> *id_mapper) {
    using namespace debruijn_graph;
    if (utils::ends_with(filename, ".gfa")) {
        gfa::GFAReader gfa(filename);
        INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links());
        gfa.to_graph(graph, id_mapper);
    } else {
        io::binary::Load(filename, graph);
    }
}

int main(int argc, char* argv[]) {
    utils::segfault_handler sh;
    gcfg cfg;

    srand(42);
    srandom(42);

    process_cmdline(argc, argv, cfg);

    create_console_logger();

    START_BANNER("SPAdes sequence-to-graph mapper");

    cfg.nthreads = spades_set_omp_threads(cfg.nthreads);
    INFO("Maximum # of threads to use (adjusted due to OMP capabilities): " << cfg.nthreads);

    try {
        unsigned nthreads = cfg.nthreads;
        unsigned k = cfg.k;
        std::string tmpdir = cfg.tmpdir, dataset_desc = cfg.file;

        fs::make_dir(tmpdir);

        DataSet dataset;
        dataset.load(dataset_desc);

        CHECK_FATAL_ERROR(cfg.libindex < dataset.lib_count(), "invalid library index");

        debruijn_graph::GraphPack gp(k, tmpdir, dataset.lib_count());
        std::unique_ptr<io::IdMapper<std::string>> id_mapper(new io::IdMapper<std::string>());

        const auto &graph = gp.get<Graph>();
        INFO("Loading de Bruijn graph from " << cfg.graph);
        LoadGraph(gp.get_mutable<Graph>(), cfg.graph, id_mapper.get());
        INFO("Graph loaded. Total vertices: " << graph.size() << ", total edges: " << graph.e_size());

        // FIXME: Get rid of this "/" junk
        debruijn_graph::config::init_libs(dataset, nthreads, tmpdir + "/");

        auto &lib = dataset[cfg.libindex];

        if (lib.is_contig_lib() || lib.is_long_read_lib() || lib.is_single()) {
            auto& path_storage = gp.get_mutable<LongReadContainer<Graph>>()[cfg.libindex];
            auto& trusted_paths = gp.get_mutable<path_extend::TrustedPathsContainer>()[cfg.libindex];

            if (lib.is_contig_lib() || lib.is_single()) {
                INFO("Mapping sequencing library #" << cfg.libindex);
                ProcessContigs(graph, cfg.mode,
                               lib,
                               path_storage,
                               trusted_paths);
            } else {
                gap_closing::GapStorage gap_storage(graph);
                PacbioAlignLibrary(graph, lib,
                                   path_storage, gap_storage,
                                   nthreads, debruijn_graph::config::pacbio_processor());
            }

            INFO("Saving to " << cfg.outfile);

            std::ofstream os(cfg.outfile);
            //FIXME fix behavior when we don't have the mapper
            path_extend::GFAPathWriter gfa_writer(graph, os,
                                                  io::MapNamingF<debruijn_graph::ConjugateDeBruijnGraph>(*id_mapper));
            gfa_writer.WriteSegmentsAndLinks();

            std::vector<PathInfo<Graph>> paths;
            path_storage.SaveAllPaths(paths);
            size_t idx = 0;
            for (const auto& entry : paths) {
                idx += 1;
                gfa_writer.WritePaths(entry.path(), std::string("PATH_") + std::to_string(idx) + "_length_" + std::to_string(entry.path().size()) + "_weigth_" + std::to_string(entry.weight()),
                                      "Z:W:" + std::to_string(entry.weight()));
            }
        } else if (lib.is_paired()) {
            paired_info::PairedIndex index(graph);
            alignment::BWAReadMapper<Graph> mapper(graph, cfg.mode);

            std::unique_ptr<ThreadPool::ThreadPool> pool;
            if (cfg.nthreads > 1)
                pool = std::make_unique<ThreadPool::ThreadPool>(cfg.nthreads);
            io::ReadConverter::ConvertToBinary(lib, pool.get());
            paired_info::FillPairedIndex(graph,
                                         mapper,
                                         lib, index, { }, 0, std::numeric_limits<unsigned>::max());

            INFO("Saving to " << cfg.outfile);

            std::ofstream os(cfg.outfile);
            for (EdgeId e : graph.edges()) {
                for (auto entry : index.Get(e)) {
                    for (const auto &point : entry.second)
                        os << (*id_mapper)[e.int_id()] << "\t" << (*id_mapper)[entry.first.int_id()] << "\t" << point << "\n";
                }
            }
        }
    } catch (const std::string &s) {
        std::cerr << s << std::endl;
        return EINTR;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EINTR;
    }

    return 0;
}
