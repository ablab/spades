//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "alignment/bwa_index.hpp"
#include "alignment/bwa_sequence_mapper.hpp"
#include "alignment/long_read_mapper.hpp"
#include "alignment/sequence_mapper_notifier.hpp"
#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"
#include "configs/aligner_config.hpp"
#include "io/binary/graph.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "io/dataset_support/read_converter.hpp"
#include "io/graph/gfa_reader.hpp"
#include "paired_info/index_point.hpp"
#include "paired_info/paired_info_utils.hpp"
#include "projects/spades/hybrid_aligning.hpp"
#include "projects/spades/hybrid_gap_closer.hpp"
#include "threadpool/threadpool.hpp"
#include "utils/filesystem/temporary.hpp"
#include "utils/logger/log_writers.hpp"
#include "utils/logger/logger.hpp"
#include "utils/segfault_handler.hpp"

#include "version.hpp"

#include <clipp/clipp.h>

#include <limits>
#include <string>

void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

using namespace debruijn_graph;

struct gcfg {
    gcfg()
        : k(-1U),
          nthreads(omp_get_max_threads() / 2 + 1),
          libindex(0),
          mode(alignment::BWAIndex::AlignmentMode::Default),
          retain(alignment::BWAIndex::RetainAlignments::Default),
          hic(false)
    {}

    unsigned k;
    std::filesystem::path file;
    std::filesystem::path graph;
    std::filesystem::path tmpdir;
    std::filesystem::path outfile;
    unsigned nthreads;
    unsigned libindex;
    alignment::BWAIndex::AlignmentMode mode;
    alignment::BWAIndex::RetainAlignments retain;
    bool hic;
    bool bin_load = false;
};

void process_cmdline(int argc, char **argv, gcfg &cfg) {
  using namespace clipp;

  std::string file;
  std::string graph;
  std::string tmpdir;
  std::string outfile;

  auto cli = (
      file << value("dataset description (in YAML)"),
      graph << value("graph (in GFA)"),
      outfile << value("output filename"),
      (option("-l") & integer("value", cfg.libindex)) % "library index (0-based, default: 0)",
      (option("-k") & integer("value", cfg.k)) % "k-mer length to use",
      (option("-t") & integer("value", cfg.nthreads)) % "# of threads to use",
      (option("--tmp-dir") & value("dir", tmpdir)) % "scratch directory to use",
      (option("--bin-load").set(cfg.bin_load)) % "load binary-converted reads from tmpdir (developer option",
      (option("--hic").set(cfg.hic)) % "enable HiC-aware paired-end processing (implies -Xhic unless mode is specified)",
      (with_prefix("-X",
                   option("illumina").set(cfg.mode, alignment::BWAIndex::AlignmentMode::Illumina) |
                   option("pacbio").set(cfg.mode, alignment::BWAIndex::AlignmentMode::PacBio) |
                   option("intractg").set(cfg.mode, alignment::BWAIndex::AlignmentMode::IntraCtg) |
                   option("hic").set(cfg.mode, alignment::BWAIndex::AlignmentMode::HiC)) % "inner alignment mode (for paired-end / contigs)"),
      (with_prefix("-R",
                   option("default").set(cfg.retain, alignment::BWAIndex::RetainAlignments::Default) |
                   option("all").set(cfg.retain, alignment::BWAIndex::RetainAlignments::All) |
                   option("primary").set(cfg.retain, alignment::BWAIndex::RetainAlignments::OnlyPrimary) |
                   option("quality").set(cfg.retain, alignment::BWAIndex::RetainAlignments::QualityPrimary)) % "retain alignments (for paired-end / contigs)")
  );

  auto result = parse(argc, argv, cli);
  if (!result) {
      auto fmt = doc_formatting{}
                 .merge_alternative_flags_with_common_prefix(true);

      std::cout << make_man_page(cli, argv[0], fmt);
      exit(1);
  }

  cfg.file = file;
  cfg.graph = graph;
  cfg.tmpdir = tmpdir.empty() ? "tmp" : tmpdir;
  cfg.outfile = outfile.empty() ? "-" : outfile;

  if (cfg.hic) {
      if (cfg.mode == alignment::BWAIndex::AlignmentMode::Default)
          cfg.mode = alignment::BWAIndex::AlignmentMode::HiC;
      
      if (cfg.retain == alignment::BWAIndex::RetainAlignments::Default)
          cfg.retain = alignment::BWAIndex::RetainAlignments::QualityPrimary;
  }
}

typedef io::DataSet<debruijn_graph::config::LibraryData> DataSet;
typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> SequencingLib;

static void ProcessContigs(const Graph &graph,
                           alignment::BWAIndex::AlignmentMode mode, alignment::BWAIndex::RetainAlignments retain,
                           SequencingLib &lib,
                           PathStorage<Graph> &single_long_reads,
                           path_extend::GappedPathStorage &trusted_paths) {
    SequenceMapperNotifier notifier;
    LongReadMapper read_mapper(graph, single_long_reads, trusted_paths, lib.type());
    notifier.Subscribe(&read_mapper);

    INFO("Mapping using BWA-mem mapper");
    alignment::BWAReadMapper<Graph> mapper(graph, mode, retain);
    auto single_streams = single_easy_readers(lib, false, false, /*handle Ns*/ false);
    notifier.ProcessLibrary(single_streams, mapper);
}

void LoadGraph(debruijn_graph::ConjugateDeBruijnGraph &graph, const std::filesystem::path &filename,
               io::IdMapper<std::string> *id_mapper) {
    using namespace debruijn_graph;
    if (filename.extension() == ".gfa") {
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
        std::filesystem::path tmpdir = cfg.tmpdir, dataset_desc = cfg.file;

        create_directory(tmpdir);

        DataSet dataset;
        dataset.load(dataset_desc);

        CHECK_FATAL_ERROR(cfg.libindex < dataset.lib_count(), "invalid library index");

        std::unique_ptr<io::IdMapper<std::string>> id_mapper(new io::IdMapper<std::string>());
        std::unique_ptr<gfa::GFAReader> gfa;
        INFO("Loading de Bruijn graph from " << cfg.graph);
        if (utils::ends_with(cfg.graph, ".gfa")) {
            gfa.reset(new gfa::GFAReader(cfg.graph));
            INFO("GFA segments: " << gfa->num_edges() << ", links: " << gfa->num_links() << ", paths: " << gfa->num_paths());
            INFO("Detected k:" << gfa->k());
            VERIFY_MSG(gfa->k() != -1U, "Failed to determine k-mer length");
            VERIFY_MSG(gfa->k() == 0 || gfa->k() % 2 == 1, "k-mer length must be odd");
            k = gfa->k();
        } else if (cfg.k == -1U)
            FATAL_ERROR("k-mer length should be specified");

        Graph graph(k);
        if (gfa) {
            gfa->to_graph(graph, id_mapper.get());
        } else {
            io::binary::Load(cfg.graph, graph);
        }
        INFO("Graph loaded. Total vertices: " << graph.size() << ", total edges: " << graph.e_size());

        debruijn_graph::config::init_libs(dataset, nthreads, tmpdir);

        auto &lib = dataset[cfg.libindex];

        if (lib.is_contig_lib() || lib.is_long_read_lib() || lib.is_single()) {
            PathStorage<Graph> path_storage(graph);

            if (lib.is_contig_lib() || lib.is_single()) {
                path_extend::GappedPathStorage trusted_paths;

                INFO("Mapping sequencing library #" << cfg.libindex);
                ProcessContigs(graph, cfg.mode, cfg.retain,
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
                gfa_writer.WritePaths(entry.path(),
                                      std::string("PATH_") + std::to_string(idx) + "_length_" + std::to_string(entry.path().size()) + "_weigth_" + std::to_string(entry.weight()),
                                      "Z:W:" + std::to_string(entry.weight()));
            }
        } else if (lib.is_paired()) {
            paired_info::PairedIndex index(graph);
            alignment::BWAReadMapper<Graph> mapper(graph, cfg.mode, cfg.retain);

            std::unique_ptr<ThreadPool::ThreadPool> pool;
            if (cfg.nthreads > 1)
                pool = std::make_unique<ThreadPool::ThreadPool>(cfg.nthreads);
            if (!cfg.bin_load || !io::ReadConverter::LoadLibIfExists(lib))
                io::ReadConverter::ConvertToBinary(lib, pool.get());

            paired_info::FillPairedIndex(graph,
                                         mapper,
                                         lib, index, { }, 0, std::numeric_limits<unsigned>::max());

            INFO("Saving to " << cfg.outfile);

            std::ofstream os(cfg.outfile);
            if (cfg.hic) {
                for (EdgeId e1 : graph.canonical_edges()) {
                    for (auto entry : index.GetHalf(e1)) {
                        EdgeId e2 = entry.first, ce2 = graph.conjugate(e2);
                        VERIFY(entry.second.size() == 1);
                        if (!(e2 <= ce2))
                            continue;

                        omnigraph::de::DEWeight w = 0;
                        auto AddToWeight = [&](EdgeId e1, EdgeId e2) {
                            const auto& hist = index.Get(e1, e2);
                            if (!hist.empty())
                                w += hist.begin()->weight;
                        };


                        if (e2 != e1) {
                            // Need to aggregate links:
                            // e1 => e2, e1 => e2', e2' => e1, e2 => e1
                            AddToWeight(e1, e2); AddToWeight(e1, ce2);
                            AddToWeight(e2, e1); AddToWeight(ce2, e1);
                        } else { // e2 == e1
                            // Need to aggregate links:
                            // e1 => e1, e1' => e1, e1 => e1'
                            AddToWeight(e1, e1); AddToWeight(ce2, e1); AddToWeight(e1, ce2);
                        }

                        os << (*id_mapper)[e1.int_id()] << "\t" << e1 << "\t" << (*id_mapper)[e2.int_id()] << "\t" << e2 << "\t" << w << "\n";
                    }
                }
            } else {
                for (EdgeId e : graph.edges()) {
                    for (auto entry : index.GetHalf(e)) {
                        for (const auto &point : entry.second)
                            os << (*id_mapper)[e.int_id()] << "\t" << (*id_mapper)[entry.first.int_id()] << "\t" << point << "\n";
                    }
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
