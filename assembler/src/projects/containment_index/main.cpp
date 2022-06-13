//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "barcode_index_construction.hpp"
#include "long_edge_statistics.hpp"
#include "multiplex_gfa_reader.hpp"
#include "reference_path_checker.hpp"
#include "scaffold_graph_helper.hpp"
#include "vertex_info_getter.hpp"
#include "vertex_resolver.hpp"
#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"

#include "io/binary/read_cloud.hpp"
#include "io/graph/gfa_writer.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/edge_cluster_extractor.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/initial_cluster_storage_builder.hpp"
#include "modules/path_extend/read_cloud_path_extend/fragment_statistics/distribution_extractor_helper.hpp"
#include "modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_polisher.hpp"
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

struct gcfg {
  size_t k = 55;
  unsigned mapping_k = 31;
  std::filesystem::path graph;
  std::filesystem::path output_dir;
  std::filesystem::path refpath;
  unsigned nthreads = (omp_get_max_threads() / 2 + 1);
  std::filesystem::path file = "";
  std::filesystem::path tmpdir = "saves";
  unsigned libindex = -1u;
  double score_threshold = 3.0;
  bool bin_load = false;
  bool debug = false;
  GraphType graph_type = GraphType::Blunted;
};

static void process_cmdline(int argc, char** argv, gcfg& cfg) {
    using namespace clipp;

    std::string graph;
    std::string output_dir;
    std::string refpath;
    std::string file;
    std::string tmpdir;

    auto cli = (
        graph << value("graph (in binary or GFA)"),
        output_dir << value("path to output directory"),
        (option("--dataset") & value("yaml", file)) % "dataset description (in YAML)",
        (option("-l") & integer("value", cfg.libindex)) % "library index (0-based, default: 0)",
        (option("-t") & integer("value", cfg.nthreads)) % "# of threads to use",
        (option("--mapping-k") & integer("value", cfg.mapping_k)) % "k for read mapping",
        (option("--tmp-dir") & value("tmp", tmpdir)) % "scratch directory to use",
        (option("--bin-load").set(cfg.bin_load)) % "load binary-converted reads from tmpdir (developer option)",
        (option("--debug").set(cfg.debug)) % "produce lots of debug data (developer option)",
        (option("--score") & value("score", cfg.score_threshold)) % "Score threshold for link index",
        (option("--ref") & value("reference", refpath)) % "Reference path",
        (with_prefix("-G",
                     option("lja").set(cfg.graph_type, GraphType::Multiplexed) |
                     option("blunt").set(cfg.graph_type, GraphType::Blunted)) % "assembly graph type")
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

void GetLongEdgeStatistics(const Graph &graph,
                           const barcode_index::FrameBarcodeIndex<Graph> &barcode_index,
                           size_t training_edge_length,
                           size_t training_edge_offset,
                           size_t read_count_threshold,
                           size_t read_linkage_distance,
                           size_t max_threads,
                           const std::filesystem::path &base_output_path) {
    cont_index::LongEdgeStatisticsCounter long_edge_counter(graph, barcode_index, training_edge_length,
                                                            training_edge_offset,
                                                            read_count_threshold, read_linkage_distance, max_threads,
                                                            base_output_path);
    long_edge_counter.CountClusterStatistics();
    long_edge_counter.CountDoubleCoverageDistribution();
}

void GetPathClusters(const Graph &graph,
                     std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr,
                     const scaffold_graph::ScaffoldGraph &scaffold_graph,
                     size_t read_linkage_distance,
                     double relative_score_threshold,
                     size_t min_read_threshold,
                     size_t length_threshold,
                     size_t max_threads,
                     io::IdMapper<std::string> *id_mapper,
                     bool lja,
                     const std::filesystem::path &output_dir) {
    INFO("Constructing initial cluster storage");
    path_extend::read_cloud::PathExtractionParams path_extraction_params(read_linkage_distance,
                                                                         relative_score_threshold,
                                                                         min_read_threshold,
                                                                         length_threshold);
    std::set<scaffold_graph::ScaffoldVertex> target_edges;
    std::copy(scaffold_graph.vbegin(), scaffold_graph.vend(), std::inserter(target_edges, target_edges.begin()));
    auto edge_cluster_extractor =
        std::make_shared<cluster_storage::AccurateEdgeClusterExtractor>(graph, barcode_extractor_ptr,
                                                                        read_linkage_distance, min_read_threshold);
    auto storage_builder =
        std::make_shared<cluster_storage::EdgeInitialClusterStorageBuilder>(graph, edge_cluster_extractor,
                                                                            target_edges, read_linkage_distance,
                                                                            min_read_threshold, max_threads);
    auto storage =
        std::make_shared<cluster_storage::InitialClusterStorage>(storage_builder->ConstructInitialClusterStorage());
    path_extend::read_cloud::ScaffoldGraphPathClusterHelper path_extractor_helper(graph,
                                                                                  barcode_extractor_ptr,
                                                                                  storage,
                                                                                  read_linkage_distance,
                                                                                  max_threads);
    INFO(scaffold_graph.VertexCount() << " vertices and " << scaffold_graph.EdgeCount() << " edges in scaffold graph");

    INFO(storage->get_cluster_storage().Size() << " initial clusters");
    auto all_clusters = path_extractor_helper.GetAllClusters(scaffold_graph);
    INFO(all_clusters.size() << " total clusters");
    std::map<size_t, size_t> size_hist;
    for (const auto &cluster: all_clusters) {
        size_hist[cluster.Size()]++;
    }
    std::ofstream sizes(output_dir / "cluster_sizes.tsv");
    for (const auto &entry: size_hist) {
        sizes << entry.first << "\t" << entry.second << std::endl;
    }
    auto path_clusters = path_extractor_helper.GetPathClusters(all_clusters);
    INFO(path_clusters.size() << " path clusters");
    size_hist.clear();
    for (const auto &cluster: path_clusters) {
        size_hist[cluster.Size()]++;
    }
    std::ofstream path_sizes(output_dir / "path_sizes.tsv");
    for (const auto &entry: size_hist) {
        path_sizes << entry.first << "\t" << entry.second << std::endl;
    }

    //extract paths traversing complex vertices
    if (lja) {
        std::unordered_set<scaffold_graph::ScaffoldVertex> uncompressed_vertices;
        for (const auto &vertex: scaffold_graph.vertices()) {
            if (not id_mapper->count(vertex.int_id())) {
                uncompressed_vertices.insert(vertex);
            }
        }

//        std::ofstream link_stream(output_dir / "lja_links.tsv");
//        for (const auto &entry: link_to_weight) {
//            link_stream << entry.first.first << "\t" << entry.first.second << "\t" << entry.second << std::endl;
//        }
    }

}

struct ScoreEntry {
  ScoreEntry(const scaffold_graph::ScaffoldVertex &first,
             const scaffold_graph::ScaffoldVertex &second,
             double hifi_score,
             double tellseq_score) : first_(first),
                                     second_(second),
                                     hifi_score_(hifi_score),
                                     tellseq_score_(tellseq_score) {}

  scaffold_graph::ScaffoldVertex first_;
  scaffold_graph::ScaffoldVertex second_;
  double hifi_score_;
  double tellseq_score_;
};

void CompareLinks(const scaffold_graph::ScaffoldGraph &hifi_graph,
                  const scaffold_graph::ScaffoldGraph &tellseq_graph,
                  cont_index::LinkIndexGraphConstructor::BarcodeScoreFunctionPtr score_function,
                  io::IdMapper<std::string> *id_mapper,
                  const std::string &output_path) {

    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    std::set<std::pair<ScaffoldVertex, ScaffoldVertex>> hifi_links;
    std::map<std::pair<ScaffoldVertex, ScaffoldVertex>, double> hifi_score_map;
    double score_threshold = 5.0;
    INFO(hifi_graph.EdgeCount() << " edges in hifi graph");
    INFO(tellseq_graph.EdgeCount() << " edges in tellseq graph");
    for (const auto &edge: hifi_graph.edges()) {
        std::pair<ScaffoldVertex, ScaffoldVertex> vertex_pair(edge.getStart(), edge.getEnd());
        hifi_links.emplace(edge.getStart(), edge.getEnd());
        hifi_score_map.emplace(vertex_pair, edge.getWeight());
    }

    std::vector<ScoreEntry> score_entries;

    if (tellseq_graph.VertexCount() != 0) {
        size_t new_tellseq_links = 0;
        for (const auto &edge: tellseq_graph.edges()) {
            if (math::ge(edge.getWeight(), score_threshold)) {
                std::pair<ScaffoldVertex, ScaffoldVertex> link(edge.getStart(), edge.getEnd());
                if (hifi_links.find(link) == hifi_links.end()) {
                    new_tellseq_links++;
                } else {
                    score_entries.emplace_back(edge.getStart(), edge.getEnd(), hifi_score_map.at(link), edge.getWeight());
                }
            }
        }
        INFO("Found " << new_tellseq_links << " new tellseq links");
    } else {
        for (const auto &edge: hifi_graph.edges()) {
            double tellseq_score = score_function->GetScore(edge);
            if (math::ge(tellseq_score, score_threshold)) {
                score_entries.emplace_back(edge.getStart(), edge.getEnd(), edge.getWeight(), tellseq_score);
            }
        }
    }
    INFO(score_entries.size() << " common links in hifi and tellseq graphs");
    std::ofstream os(output_path);
    os << "First\tSecond\tHiFi links\tTellSeq links" << "\n";
    for (const auto &entry: score_entries) {
        os << (*id_mapper)[entry.first_.int_id()] << "\t" << (*id_mapper)[entry.second_.int_id()]
           << "\t" << entry.hifi_score_ << "\t" << entry.tellseq_score_ << "\n";
    }
}

void NormalizeTellseqLinks(const scaffold_graph::ScaffoldGraph &tellseq_graph,
                           size_t min_length,
                           std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr,
                           size_t count_threshold,
                           size_t tail_threshold,
                           io::IdMapper<std::string> *id_mapper,
                           const std::filesystem::path &output_path) {
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    const auto &assembly_graph = tellseq_graph.AssemblyGraph();
    std::map<std::pair<ScaffoldVertex, ScaffoldVertex>, double> score_map;
    double score_threshold = 2.0;
    if (tellseq_graph.VertexCount() != 0) {
        for (const auto &edge: tellseq_graph.edges()) {
            const auto &first_vertex = edge.getStart();
            const auto &first_conj = first_vertex.GetConjugateFromGraph(assembly_graph);
            const auto &second_vertex = edge.getEnd();
            const auto &second_conj = second_vertex.GetConjugateFromGraph(assembly_graph);
            size_t first_len = first_vertex.GetLengthFromGraph(assembly_graph);
            size_t second_len = second_vertex.GetLengthFromGraph(assembly_graph);
            double num_barcodes = edge.getWeight();
            if (first_len >= min_length and second_len >= min_length and math::ge(num_barcodes, score_threshold)) {
                std::pair<ScaffoldVertex, ScaffoldVertex> link(first_vertex, second_vertex);
                std::pair<ScaffoldVertex, ScaffoldVertex> rc_rc_link(second_conj, first_conj);
                if (first_vertex != second_vertex and first_vertex != second_conj) {
                    score_map[link] += num_barcodes / 2;
                    score_map[rc_rc_link] += num_barcodes / 2;
                }
            }
        }
    }
    INFO(score_map.size() << " filtered tellseq links");

    std::unordered_map<ScaffoldVertex, size_t> vertex_to_head_barcodes;
    for (const auto &vertex: tellseq_graph.vertices()) {
        vertex_to_head_barcodes[vertex] = barcode_extractor_ptr->GetBarcodesFromHead(vertex.GetFirstEdge(), count_threshold, tail_threshold).size();
    }

    std::ofstream os(output_path / "normalized_tellseq_links.tsv");
    os << "First\tSecond\tTotal links\tJaccard Index" << "\n";
    for (const auto &entry: score_map) {
        const auto &first_vertex = entry.first.first;
        const auto &first_conj = first_vertex.GetConjugateFromGraph(assembly_graph);
        const auto &second_vertex = entry.first.second;
        double shared_barcodes = entry.second;
        auto first_tail_barcodes = vertex_to_head_barcodes[first_conj];
        auto second_head_barcodes = vertex_to_head_barcodes[second_vertex];
        auto union_size = static_cast<double>(first_tail_barcodes + second_head_barcodes) - shared_barcodes;
        double jaccard_index = shared_barcodes / union_size;
        os << (*id_mapper)[first_vertex.int_id()] << "\t" << (*id_mapper)[second_vertex.int_id()]
           << "\t" << shared_barcodes << "\t" << jaccard_index << "\n";
    }
}

void ConstructCloudOnlyLinks(const Graph &graph,
                             std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr,
                             double graph_score_threshold,
                             size_t length_threshold,
                             size_t tail_threshold,
                             size_t count_threshold,
                             size_t threads,
                             bool bin_load,
                             bool debug,
                             const std::filesystem::path &output_dir,
                             io::IdMapper<std::string> *id_mapper) {
    auto tellseq_graph = cont_index::GetTellSeqScaffoldGraph(graph, barcode_extractor_ptr, graph_score_threshold,
                                                             length_threshold, tail_threshold, count_threshold,
                                                             threads, bin_load, debug, output_dir, id_mapper);

//    LinkIndexGraphConstructor link_index_constructor(graph, barcode_extractor_ptr, graph_score_threshold,
//                                                     tail_threshold, length_threshold, count_threshold, threads);
//    auto score_function = link_index_constructor.ConstructScoreFunction();
//    NormalizeTellseqLinks(tellseq_graph, length_threshold, barcode_extractor_ptr, count_threshold,
//                          tail_threshold, id_mapper, output_dir);
//    auto compare_output_path = output_dir / "hifi_tellseq_scores.tsv";
//    CompareLinks(hifi_graph, tellseq_graph, score_function, id_mapper.get(), compare_output_path);

}


void ReadGraph(const gcfg &cfg,
               debruijn_graph::Graph &graph,
               io::IdMapper<std::string> *id_mapper) {
    switch (cfg.graph_type) {
        default:
            FATAL_ERROR("Unknown graph representation type");
        case GraphType::Blunted: {
            gfa::GFAReader gfa(cfg.graph);
            INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links() << ", paths: "
                                  << gfa.num_paths());
            gfa.to_graph(graph, id_mapper);
            return;
        }
        case GraphType::Multiplexed: {
//            cont_index::MultiplexGFAReader demulti_gfa(cfg.graph);
//            gfa::GFAReader gfa(cfg.graph);
//            INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links() << ", paths: " << gfa.num_paths());
//            VERIFY_MSG(demulti_gfa.k() != -1U, "Failed to determine k-mer length");
//            INFO(demulti_gfa.k());
//            debruijn_graph::Graph struct_graph(demulti_gfa.k());
//            gfa.to_graph(struct_graph, struct_id_mapper.get());
//            debruijn_graph::Graph graph(demulti_gfa.k());
//            demulti_gfa.to_graph(struct_graph, struct_id_mapper.get(), graph, id_mapper.get());
            FATAL_ERROR("Multiplexed graphs are currently not supported");
        }
    }
}

void AnalyzeVertices(debruijn_graph::Graph &graph,
                     std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr,
                     size_t count_threshold,
                     size_t tail_threshold,
                     size_t length_threshold,
                     size_t threads,
                     double score_threshold,
                     io::IdMapper<std::string> *id_mapper,
                     const std::filesystem::path &output_path,
                     const std::filesystem::path &tmp_path) {
    cont_index::VertexResolver vertex_resolver
        (graph, barcode_extractor_ptr, count_threshold, tail_threshold, length_threshold, threads, score_threshold);
    const auto &vertex_results = vertex_resolver.ResolveVertices() ;
    std::filesystem::path vertex_output_path = output_path / "vertex_stats.tsv";
    vertex_resolver.PrintVertexResults(vertex_results, vertex_output_path, tmp_path, id_mapper);
    cont_index::PathExtractor path_extractor(graph);
    const auto &paths = path_extractor.ExtractPaths(vertex_results);
    path_extractor.PrintPaths(paths, output_path / "contigs.paths", id_mapper);
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
    gfa::GFAWriter gfa_writer(graph, graph_out);
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

        //fixme configs
        //index construction
        const size_t frame_size = 2000;
        const size_t read_linkage_distance = 40000;

        //graph construction
        const double graph_score_threshold = 1.99;
        const size_t tail_threshold = 200000;
        const size_t length_threshold = 0;
        const size_t count_threshold = 1;

        //path cluster extraction
        const double relative_score_threshold = 10.0;
        const size_t min_read_threshold = 2;

        //Long edge clusters
        const size_t training_length_threshold = 500000;
        const size_t training_length_offset = 50000;
        const size_t training_min_read_threshold = 6;
        const size_t training_linkage_distance = 30000;

        debruijn_graph::config::init_libs(dataset, cfg.nthreads, cfg.tmpdir);
        barcode_index::FrameBarcodeIndex<debruijn_graph::Graph> barcode_index(graph, read_linkage_distance);
        using BarcodeExtractor = barcode_index::FrameBarcodeIndexInfoExtractor;
        auto barcode_extractor_ptr = std::make_shared<BarcodeExtractor>(barcode_index, graph);

        std::unique_ptr<TimeTracerRAII> traceraii;
        traceraii.reset(new TimeTracerRAII(argv[0], 500));
        INFO("Time tracing is enabled");

        TIME_TRACE_SCOPE("Containment index");

        auto &lib = dataset[cfg.libindex];
        if (lib.type() == io::LibraryType::Clouds10x) {
            cont_index::ConstructBarcodeIndex(barcode_index,
                                              lib,
                                              graph,
                                              cfg.tmpdir,
                                              cfg.nthreads,
                                              frame_size,
                                              cfg.mapping_k,
                                              cfg.bin_load,
                                              cfg.debug);
        } else {
            WARN("Only read cloud libraries with barcode tags are supported for links");
        }

        AnalyzeVertices(graph, barcode_extractor_ptr, count_threshold, tail_threshold, length_threshold, cfg.nthreads,
                        2.0, id_mapper.get(), cfg.output_dir, cfg.tmpdir);


//        GetLongEdgeStatistics(graph, barcode_index, training_length_threshold, training_length_offset,
//                              training_min_read_threshold, training_linkage_distance, cfg.nthreads, cfg.output_dir);

//        if (not cfg.refpath.empty()) {
//            cont_index::ReferencePathChecker ref_path_checker(graph, id_mapper.get());
//            ref_path_checker.CheckAssemblyGraph(cfg.refpath);
//        }
        bool lja = true;
//        GetPathClusters(graph, barcode_extractor_ptr, hifi_graph, read_linkage_distance, relative_score_threshold,
//                        min_read_threshold, length_threshold, cfg.nthreads, id_mapper.get(), lja, cfg.output_dir);
//        ConstructCloudOnlyLinks(graph, barcode_extractor_ptr, graph_score_threshold, length_threshold, tail_threshold,
//                                count_threshold, cfg.nthreads, cfg.bin_load, cfg.debug, cfg.output_dir, id_mapper.get());
    }
}
