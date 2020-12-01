//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"

#include "modules/alignment/sequence_mapper_notifier.hpp"
#include "modules/alignment/bwa_sequence_mapper.hpp"
#include "modules/alignment/long_read_mapper.hpp"

#include "io/graph/gfa_reader.hpp"
#include "io/dataset_support/read_converter.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "io/binary/graph.hpp"

#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"

#include "utils/logger/log_writers.hpp"
#include "utils/logger/logger.hpp"
#include "utils/segfault_handler.hpp"
#include "utils/filesystem/temporary.hpp"

#include "pipeline/graph_pack.hpp"
#include "pipeline/configs/aligner_config.hpp"

#include "projects/spades/hybrid_aligning.hpp"
#include "projects/spades/hybrid_gap_closer.hpp"

#include "version.hpp"

#include <clipp/clipp.h>

#include <sys/types.h>
#include <sys/stat.h>

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
      (option("--tmp-dir") & value("dir", cfg.tmpdir)) % "scratch directory to use"
  );

  auto result = parse(argc, argv, cli);
  if (!result) {
      std::cout << make_man_page(cli, argv[0]);
      exit(1);
  }
}

typedef io::DataSet<debruijn_graph::config::LibraryData> DataSet;
typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> SequencingLib;

std::shared_ptr<SequenceMapper<Graph>> ChooseProperMapper(const GraphPack& gp,
                                                          const SequencingLib& library) 
{
    const auto &graph = gp.get<Graph>();
    if (library.type() == io::LibraryType::MatePairs) {
        INFO("Mapping mate pairs using BWA-mem mapper");
        return std::make_shared<alignment::BWAReadMapper<Graph>>(graph);
    }

    if (library.data().unmerged_read_length < gp.k() && library.type() == io::LibraryType::PairedEnd) {
        INFO("Mapping PE reads shorter than K with BWA-mem mapper");
        return std::make_shared<alignment::BWAReadMapper<Graph>>(graph);
    }

    INFO("Selecting usual mapper");
    return MapperInstance(gp);
}

static void ProcessSingleReads(GraphPack &gp, DataSet &dataset,
                               size_t ilib,
                               bool use_binary = true,
                               bool map_paired = false) {
    SequencingLib &lib = dataset[ilib];
    SequenceMapperNotifier notifier(gp, dataset.lib_count());

    const auto &graph = gp.get<Graph>();
    auto& single_long_reads = gp.get_mutable<LongReadContainer<Graph>>()[ilib];
    auto& trusted_paths = gp.get_mutable<path_extend::TrustedPathsContainer>()[ilib];
    LongReadMapper read_mapper(graph, single_long_reads, trusted_paths, lib.type());

    if (lib.is_contig_lib()) {
        //FIXME pretty awful, would be much better if listeners were shared ptrs
        notifier.Subscribe(ilib, &read_mapper);
    }

    auto mapper_ptr = ChooseProperMapper(gp, lib);
    if (use_binary) {
        auto single_streams = single_binary_readers(lib, false, map_paired);
        notifier.ProcessLibrary(single_streams, ilib, *mapper_ptr);
    } else {
        auto single_streams = single_easy_readers(lib, false,
                                                  map_paired, /*handle Ns*/false);
        notifier.ProcessLibrary(single_streams, ilib, *mapper_ptr);
    }
}

#if 0
static void ProcessPairedReads(GraphPack &gp,
                               DataSet &dataset, size_t ilib) {
    SequencingLib &lib = dataset[ilib];
    SequenceMapperNotifier notifier(gp, dataset.lib_count());

    auto paired_streams = paired_binary_readers(lib, /*followed by rc*/false, 0, /*include merged*/true);
    notifier.ProcessLibrary(paired_streams, ilib, *ChooseProperMapper(gp, lib));
}
#endif

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

        debruijn_graph::GraphPack gp(k, tmpdir, dataset.lib_count());
        std::unique_ptr<io::IdMapper<std::string>> id_mapper(new io::IdMapper<std::string>());

        const auto &graph = gp.get<Graph>();
        INFO("Loading de Bruijn graph from " << cfg.graph);
        LoadGraph(gp.get_mutable<Graph>(), cfg.graph, id_mapper.get());
        {
            size_t sz = 0;
            for (auto it = graph.ConstEdgeBegin(); !it.IsEnd(); ++it)
                sz += 1;

            INFO("Graph loaded. Total vertices: " << graph.size() << " Total edges: " << sz);
        }

        // FIXME: Get rid of this "/" junk
        debruijn_graph::config::init_libs(dataset, nthreads, tmpdir + "/");

        gp.get_mutable<KmerMapper<Graph>>().Attach();
        gp.EnsureBasicMapping();

        for (size_t i = 0; i < dataset.lib_count(); ++i) {
            auto &lib = dataset[i];
            auto& path_storage = gp.get_mutable<LongReadContainer<Graph>>()[i];
            if (lib.is_contig_lib()) {
                INFO("Mapping contigs library #" << i);
                ProcessSingleReads(gp, dataset, i, false);
            } else if (lib.is_long_read_lib()) {
                gap_closing::GapStorage gap_storage(graph);

                debruijn_graph::config::pacbio_processor pb;

                PacbioAlignLibrary(graph, lib,
                                   path_storage, gap_storage,
                                   nthreads, pb);
            } else {
                WARN("Could only map contigs or long reads so far, skipping the library");
                continue;
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
