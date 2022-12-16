//***************************************************************************
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "barcode_index_construction.hpp"
#include "cloud_only_links.hpp"
#include "graph_resolver.hpp"
#include "graph_resolver_io.hpp"
#include "long_edge_statistics.hpp"
#include "path_cluster_extractor.hpp"
#include "path_extractor.hpp"
#include "reference_path_checker.hpp"
#include "scaffold_graph_helper.hpp"
#include "vertex_info_getter.hpp"
#include "vertex_resolver.hpp"

#include "auxiliary_graphs/contracted_graph/contracted_graph_builder.hpp"
#include "io/binary/read_cloud.hpp"
#include "io/graph/gfa_writer.hpp"
#include "modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/path_cluster_helper.hpp"
#include "toolchain/utils.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "utils/parallel/openmp_wrapper.h"
#include "utils/segfault_handler.hpp"
#include "utils/verify.hpp"

#include <clipp/clipp.h>

using namespace debruijn_graph;
using namespace cont_index;
using namespace path_extend::read_cloud;

enum class GraphType {
    Blunted,
    Multiplexed
};

enum class ResolutionMode {
    Diploid,
    Meta
};

struct gcfg {
  unsigned k = 55;
  unsigned mapping_k = 31;
  std::filesystem::path graph;
  std::filesystem::path output_dir;
  std::filesystem::path refpath;
  std::filesystem::path assembly_info;
  unsigned nthreads = (omp_get_max_threads() / 2 + 1);
  std::filesystem::path file = "";
  std::filesystem::path tmpdir = "saves";
  unsigned libindex = 0;
  GraphType graph_type = GraphType::Multiplexed;
  ResolutionMode mode = ResolutionMode::Diploid;
  bool bin_load = false;
  bool debug = false;
  bool statistics = false;

  //barcode_index_construction
  size_t frame_size = 40000;
  size_t read_linkage_distance = 40000;

  //graph construction
  double graph_score_threshold = 2.0;
  size_t tail_threshold = 200000;
  size_t count_threshold = 1;

  //meta mode
  size_t length_threshold = 2000;

  //path cluster extraction
  double relative_score_threshold = 10.0;
  size_t min_read_threshold = 2;

  //Long edge clusters
  size_t training_length_threshold = 500000;
  size_t training_length_offset = 50000;
  size_t training_min_read_threshold = 6;
  size_t training_linkage_distance = 30000;
};

static void process_cmdline(int argc, char** argv, gcfg& cfg) {
    using namespace clipp;

    std::string graph;
    std::string output_dir;
    std::string refpath;
    std::string file;
    std::string tmpdir;
    std::string assembly_info;

    auto cli = (
        graph << value("graph (in binary or GFA)"),
        file << value("SLR library description (in YAML)"),
        output_dir << value("path to output directory"),
        (option("--dataset") & value("yaml", file)) % "dataset description (in YAML)",
        (option("-l") & integer("value", cfg.libindex)) % "library index (0-based, default: 0)",
        (option("-t") & integer("value", cfg.nthreads)) % "# of threads to use",
        (option("--mapping-k") & integer("value", cfg.mapping_k)) % "k for read mapping",
        (option("--tmp-dir") & value("tmp", tmpdir)) % "scratch directory to use",
        (option("--ref") & value("reference", refpath)) % "Reference path for repeat resolution evaluation (developer option)",
        (option("--assembly-info") & value("assembly-info", assembly_info))
            % "Path to metaflye assembly_info.txt file (meta mode, metaFlye graphs only)",
        (option("--bin-load").set(cfg.bin_load)) % "load binary-converted reads from tmpdir (developer option)",
        (option("--debug").set(cfg.debug)) % "produce lots of debug data (developer option)",
        (option("--statistics").set(cfg.statistics)) % "produce additional read cloud library statistics (developer option)",
        (with_prefix("-G",
                     option("mdbg").set(cfg.graph_type, GraphType::Multiplexed) |
                     option("blunt").set(cfg.graph_type, GraphType::Blunted)) % "assembly graph type (mDBG or blunted)"),
        (with_prefix("-M",
                     option("diploid").set(cfg.mode, ResolutionMode::Diploid) |
                     option("meta").set(cfg.mode, ResolutionMode::Meta)) % "repeat resolution mode (diploid or meta)"),
         (option("--frame-size") & value("frame-size", cfg.frame_size)) % "Resolution of barcode index",
        (option("--linkage-distance") & value("read-linkage-distance", cfg.read_linkage_distance)) %
            "Reads are assigned to the same fragment based on linkage distance",
        (option("--score") & value("score", cfg.graph_score_threshold)) % "Score threshold for link index",
        (option("--tail-threshold") & value("tail-threshold", cfg.tail_threshold)) %
            "Barcodes are assigned to the first and last <tail_threshold> nucleotides of the edge",
        (option("--count-threshold") & value("count-threshold", cfg.count_threshold))
            % "Minimum number of reads for barcode index",
        (option("--relative-score-threshold") & value("relative-score-threshold", cfg.relative_score_threshold))
            % "Relative score threshold for path cluster extraction",
        (option("--min-read-threshold") & value("min-read-threshold", cfg.min_read_threshold))
            % "Minimum number of reads for path cluster extraction",
        (option("--length-threshold") & value("length-threshold", cfg.length_threshold))
            % "Minimum scaffold graph edge length (meta mode option)"
    );

    auto result = parse(argc, argv, cli);
    if (!result) {
        std::cout << make_man_page(cli, argv[0]);
        exit(1);
    }
    cfg.graph = graph;
    cfg.output_dir = output_dir;
    cfg.refpath = refpath;
    cfg.file = file;
    cfg.tmpdir = tmpdir;
    cfg.assembly_info = assembly_info;
}

struct TimeTracerRAII {
  TimeTracerRAII(llvm::StringRef program_name,
                 unsigned granularity = 500,
                 const std::string &prefix = "", const std::string &suffix = "") {
      time_trace_file_ = prefix + "time_trace_" + suffix + ".json";
      llvm::timeTraceProfilerInitialize(granularity, program_name);
  }
  ~TimeTracerRAII() {
      if (auto E = llvm::timeTraceProfilerWrite(time_trace_file_, "cont-index")) {
          handleAllErrors(std::move(E),
                          [&](const llvm::StringError &SE) {
                            ERROR("" << SE.getMessage() << "\n");
                          });
          return;
      } else {
          INFO("Time trace is written to: " << time_trace_file_);
      }
      llvm::timeTraceProfilerCleanup();
  }

  std::string time_trace_file_;
};

void ReadGraph(const gcfg &cfg,
               debruijn_graph::Graph &graph,
               io::IdMapper<std::string> *id_mapper) {
    switch (cfg.graph_type) {
        default:
            FATAL_ERROR("Unknown graph representation type");
        case GraphType::Multiplexed: {
            gfa::GFAReader gfa(cfg.graph);
            INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links() << ", paths: "
                                  << gfa.num_paths());
            gfa.to_graph(graph, id_mapper);
            return;
        }
    }
}

std::unordered_set<debruijn_graph::EdgeId> ParseRepetitiveEdges(const debruijn_graph::Graph &graph,
                                                                const std::string &path_to_info,
                                                                io::IdMapper<std::string> *id_mapper) {
    std::unordered_set<EdgeId> result;
    size_t total_repetitive_length = 0;
    std::ifstream info_stream(path_to_info);
    std::string ctg_id;
    std::string is_repeat;
    std::string blank;
    std::string graph_path;
    for (size_t i = 0; i < 8; ++i) {
        info_stream >> blank;
    }
    while (!info_stream.eof()) {
        info_stream >> ctg_id;
        info_stream >> blank;
        info_stream >> blank;
        info_stream >> blank;
        info_stream >> is_repeat;
        info_stream >> blank;
        info_stream >> blank;
        info_stream >> graph_path;
//        INFO(graph_path);
        if (is_repeat == "Y") {
            auto current_pos = graph_path.find(',');
            while (current_pos != std::string::npos) {
                std::string edge_num = graph_path.substr(0, current_pos);
                if (edge_num != "*") {
                    EdgeId edge = (*id_mapper)["edge_" + edge_num];
                    result.insert(edge);
                    result.insert(graph.conjugate(edge));
                    total_repetitive_length += graph.length(edge);
                }
                graph_path.erase(0, current_pos + 1);
//                INFO(graph_path);
                current_pos = graph_path.find(',');
            }
        }
    }
    INFO(result.size() << " repetitive edges, total length: " << total_repetitive_length);
    return result;
}

size_t GetLengthThreshold(const gcfg &cfg) {
    switch (cfg.mode) {
        default:
            FATAL_ERROR("Unknown repeat resolution mode");
        case ResolutionMode::Diploid: {
            return 0;
        }
        case ResolutionMode::Meta: {
            return cfg.length_threshold;
        }
    }
}

void CountStatistics(const gcfg &cfg,
                     debruijn_graph::Graph &graph,
                     const barcode_index::FrameBarcodeIndex<Graph> &barcode_index,
                     const std::filesystem::path &base_output_path) {
    cont_index::LongEdgeStatisticsCounter long_edge_counter(graph, barcode_index, cfg.training_length_threshold,
                                                            cfg.training_length_offset,
                                                            cfg.training_min_read_threshold, cfg.read_linkage_distance,
                                                            cfg.nthreads,
                                                            base_output_path);
    long_edge_counter.CountClusterStatistics();
    long_edge_counter.CountDoubleCoverageDistribution();

    auto barcode_extractor_ptr = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(barcode_index, graph);
    cont_index::PathClusterExtractor path_cluster_extractor(barcode_extractor_ptr,
                                                            cfg.read_linkage_distance,
                                                            cfg.relative_score_threshold,
                                                            cfg.min_read_threshold,
                                                            cfg.length_threshold,
                                                            cfg.nthreads);
}

cont_index::VertexResults GetRepeatResolutionResults(const gcfg &cfg,
                                                     debruijn_graph::Graph &graph,
                                                     std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr,
                                                     io::IdMapper<std::string> *id_mapper) {
    size_t length_threshold = GetLengthThreshold(cfg);

    //fixme separate output
    std::filesystem::path vertex_output_path = cfg.output_dir / "vertex_stats.tsv";
    bool unique_kmer_count = false;
    switch (cfg.mode) {
        default:
            FATAL_ERROR("Unknown repeat resolution mode");
        case ResolutionMode::Diploid: {
            cont_index::VertexResolver<debruijn_graph::Graph> vertex_resolver
                (graph, graph, barcode_extractor_ptr, cfg.count_threshold, cfg.tail_threshold, length_threshold, cfg.nthreads,
                 cfg.graph_score_threshold);
            auto results = vertex_resolver.ResolveVertices();
            vertex_resolver.PrintVertexResults(results, vertex_output_path, cfg.tmpdir, unique_kmer_count, id_mapper);
            return results;
        }
        case ResolutionMode::Meta: {
//            auto length_predicate = [&graph, length_threshold](const debruijn_graph::EdgeId &edge) {
//              return graph.length(edge) >= length_threshold;
//            };

            auto repetitive_edges = ParseRepetitiveEdges(graph, cfg.assembly_info, id_mapper);
            auto repeat_predicate = [&repetitive_edges](const debruijn_graph::EdgeId &edge) {
                return repetitive_edges.find(edge) == repetitive_edges.end();
            };
            contracted_graph::DBGContractedGraphFactory factory(graph, repeat_predicate);
            factory.Construct();
            INFO("Constructed graph");
            auto contracted_graph = factory.GetGraph();
            INFO("Total contracted graph vertices: " << contracted_graph->size());
            cont_index::VertexResolver<contracted_graph::ContractedGraph> vertex_resolver
                (*contracted_graph, contracted_graph->GetAssemblyGraph(), barcode_extractor_ptr, cfg.count_threshold,
                 cfg.tail_threshold, length_threshold, cfg.nthreads,
                 cfg.graph_score_threshold);
            auto results = vertex_resolver.ResolveVertices();
            vertex_resolver.PrintVertexResults(results, vertex_output_path, cfg.tmpdir, unique_kmer_count, id_mapper);
            return results;
        }
    }

}

void ResolveComplexVertices(const gcfg &cfg,
                            debruijn_graph::Graph &graph,
                            std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr,
                            io::IdMapper<std::string> *id_mapper,
                            path_extend::GFAPathWriter gfa_writer) {
    auto vertex_results = GetRepeatResolutionResults(cfg, graph, barcode_extractor_ptr, id_mapper);

    cont_index::PathExtractor path_extractor(graph);
    path_extend::PathContainer paths;
    path_extractor.ExtractPaths(paths, vertex_results);
    INFO("Extracted paths")

    auto name_generator = std::make_shared<path_extend::DefaultContigNameGenerator>();
    path_extend::ContigWriter writer(graph, name_generator);
    std::vector<path_extend::PathsWriterT> path_writers;
    INFO("Creating writers")
    path_writers.push_back([&](const path_extend::ScaffoldStorage &scaffold_storage) {
      auto fn = cfg.output_dir / ("contigs.fasta");
      INFO("Outputting contigs to " << fn);
      path_extend::ContigWriter::WriteScaffolds(scaffold_storage, fn);
    });
    path_writers.push_back([&](const path_extend::ScaffoldStorage &storage) {
      INFO("Populating GFA with scaffold paths");
      gfa_writer.WritePaths(storage);
    });
    writer.OutputPaths(paths, path_writers);

    if (cfg.mode == ResolutionMode::Diploid) {
        cont_index::GraphResolver graph_resolver;
        auto graph_resolver_info = graph_resolver.TransformGraph(graph, paths, vertex_results);
        TransformedGraphIO graph_output(id_mapper);
        graph_output.PrintGraph(graph, graph_resolver_info, cfg.output_dir);
    }
}

int main(int argc, char** argv) {
    utils::segfault_handler sh;
    gcfg cfg;

    srand(42);
    srandom(42);

    process_cmdline(argc, argv, cfg);

    toolchain::create_console_logger();
    START_BANNER("Containment index builder");

    cfg.nthreads = spades_set_omp_threads(cfg.nthreads);
    INFO("Maximum # of threads to use (adjusted due to OMP capabilities): " << cfg.nthreads);

    std::filesystem::create_directory(cfg.output_dir);
    std::filesystem::create_directory(cfg.tmpdir);

    INFO("Loading graph");
    std::unique_ptr<io::IdMapper<std::string>> id_mapper(new io::IdMapper<std::string>());

    debruijn_graph::Graph graph(cfg.k);
    ReadGraph(cfg, graph, id_mapper.get());
    INFO("Graph loaded. Total vertices: " << graph.size() << ", total edges: " << graph.e_size());

    std::ofstream graph_out(cfg.output_dir / "assembly_graph.gfa");
    path_extend::GFAPathWriter gfa_writer(graph, graph_out,
                                          io::MapNamingF<debruijn_graph::ConjugateDeBruijnGraph>(*id_mapper));
    gfa_writer.WriteSegmentsAndLinks();

    INFO("Building barcode index");
    if (cfg.libindex != -1u) {
        INFO("Processing paired-end reads");
        DataSet dataset;
        if (cfg.file != "") {
            dataset.load(cfg.file);
            if (cfg.libindex == -1u)
                cfg.libindex = 0;
            CHECK_FATAL_ERROR(cfg.libindex < dataset.lib_count(), "invalid library index");
        }

        debruijn_graph::config::init_libs(dataset, cfg.nthreads, cfg.tmpdir);
        barcode_index::FrameBarcodeIndex<debruijn_graph::Graph> barcode_index(graph, cfg.frame_size);
        using BarcodeExtractor = barcode_index::FrameBarcodeIndexInfoExtractor;
        auto barcode_extractor_ptr = std::make_shared<BarcodeExtractor>(barcode_index, graph);

        std::unique_ptr<TimeTracerRAII> traceraii;
        traceraii.reset(new TimeTracerRAII(argv[0], 500));
        INFO("Time tracing is enabled");

        TIME_TRACE_SCOPE("Containment index");

        auto &lib = dataset[cfg.libindex];
        if (lib.type() == io::LibraryType::Clouds10x) {
            cont_index::ConstructBarcodeIndex(barcode_index, lib, graph, cfg.tmpdir, cfg.nthreads, cfg.frame_size,
                                              cfg.mapping_k, cfg.bin_load, cfg.debug);
        } else {
            WARN("Only read cloud libraries with barcode tags are supported for links");
        }

        ResolveComplexVertices(cfg, graph, barcode_extractor_ptr, id_mapper.get(), gfa_writer);

        if (cfg.statistics) {
            CountStatistics(cfg, graph, barcode_index, cfg.output_dir);
        }

        if (not cfg.refpath.empty()) {
            cont_index::ReferencePathChecker ref_path_checker(graph, id_mapper.get());
            ref_path_checker.CheckAssemblyGraph(cfg.refpath);
        }
    }
}
