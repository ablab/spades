//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "barcode_index_construction.hpp"
#include "multiplex_gfa_reader.hpp"
#include "scaffold_graph_helper.hpp"
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

struct gcfg {
  size_t k = 55;
  std::string graph;
  std::string output_dir;
  std::string refpath;
  unsigned nthreads = (omp_get_max_threads() / 2 + 1);
  std::string file = "";
  std::string tmpdir = "saves";
  unsigned libindex = -1u;
  double score_threshold = 3.0;
  bool bin_load = false;
  bool debug = false;
};

static void process_cmdline(int argc, char** argv, gcfg& cfg) {
    using namespace clipp;

    auto cli = (
        cfg.graph << value("graph (in binary or GFA)"),
            cfg.output_dir << value("path to output directory"),
            (option("--dataset") & value("yaml", cfg.file)) % "dataset description (in YAML)",
            (option("-l") & integer("value", cfg.libindex)) % "library index (0-based, default: 0)",
            (option("-t") & integer("value", cfg.nthreads)) % "# of threads to use",
            (option("--tmp-dir") & value("tmp", cfg.tmpdir)) % "scratch directory to use",
            (option("--bin-load").set(cfg.bin_load)) % "load binary-converted reads from tmpdir (developer option)",
            (option("--debug").set(cfg.debug)) % "produce lots of debug data (developer option)",
            (option("--score") & value("score", cfg.score_threshold)) % "Score threshold for link index",
            (option("--ref") & value("reference", cfg.refpath)) % "Reference path"
    );

    auto result = parse(argc, argv, cli);
    if (!result) {
        std::cout << make_man_page(cli, argv[0]);
        exit(1);
    }
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
                           const std::string &base_output_path) {
    auto barcode_extractor = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(barcode_index, graph);
    std::set<scaffold_graph::ScaffoldVertex> long_edges;
    for (const EdgeId &edge: graph.canonical_edges()) {
        if (graph.length(edge) >= training_edge_length) {
            long_edges.insert(edge);
        }
    }
    INFO("Found " << long_edges.size() << " edges longer than " << training_edge_length << " in the graph");
    auto edge_cluster_extractor =
        std::make_shared<cluster_storage::AccurateEdgeClusterExtractor>(graph, barcode_extractor,
                                                                        read_linkage_distance, read_count_threshold);
    cluster_storage::EdgeInitialClusterStorageBuilder initial_builder(graph, edge_cluster_extractor, long_edges,
                                                                      read_linkage_distance, read_count_threshold,
                                                                      max_threads);
    auto initial_cluster_storage = initial_builder.ConstructInitialClusterStorage();
    path_extend::read_cloud::fragment_statistics::ClusterDistributionExtractor distribution_extractor(graph, barcode_index,
                                                                                                      read_count_threshold,
                                                                                                      training_edge_length,
                                                                                                      training_edge_offset,
                                                                                                      max_threads);
    auto distribution_pack = distribution_extractor.GetDistributionsForStorage(initial_cluster_storage.get_cluster_storage());
    std::string length_output_path = fs::append_path(base_output_path, "long_edge_fragment_lengths.tsv");
    std::ofstream length_stream(length_output_path);
    length_stream << "Length\tNumber\n";
    size_t total_long_edge_clusters = 0;
    for (const auto &entry: distribution_pack.length_distribution_) {
        length_stream << entry.first << "\t" << entry.second << std::endl;
        total_long_edge_clusters += entry.second;
    }
    INFO("Total long edge clusters: " << total_long_edge_clusters);
    std::string coverage_output_path = fs::append_path(base_output_path, "long_edge_fragment_coverages.tsv");
    std::ofstream coverage_stream(coverage_output_path);
    coverage_stream << "Coverage\tNumber\n";
    for (const auto &entry: distribution_pack.coverage_distribution_) {
        coverage_stream << entry.first << "\t" << entry.second << std::endl;
    }
    std::string number_of_reads_output_path = fs::append_path(base_output_path, "long_edge_fragment_reads.tsv");
    std::ofstream read_stream(number_of_reads_output_path);
    read_stream << "Number of reads\tNumber\n";
    for (const auto &entry: distribution_pack.num_reads_distribution_) {
        read_stream << entry.first << "\t" << entry.second << std::endl;
    }
}

void GetPathClusterStatistics(const Graph &graph,
                              std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr,
                              const scaffold_graph::ScaffoldGraph &scaffold_graph,
                              size_t read_linkage_distance,
                              double relative_score_threshold,
                              size_t min_read_threshold,
                              size_t length_threshold,
                              size_t max_threads,
                              const std::string &output_dir) {
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
                                                                            min_read_threshold,
                                                                            max_threads);
    auto storage =
        std::make_shared<cluster_storage::InitialClusterStorage>(storage_builder->ConstructInitialClusterStorage());
    path_extend::read_cloud::ScaffoldGraphPathClusterHelper path_extractor_helper(graph,
                                                                                  barcode_extractor_ptr,
                                                                                  storage,
                                                                                  read_linkage_distance,
                                                                                  max_threads);

    INFO(storage->get_cluster_storage().Size() << " initial clusters");
    auto all_clusters = path_extractor_helper.GetAllClusters(scaffold_graph);
    INFO(all_clusters.size() << " total clusters");
    std::map<size_t, size_t> size_hist;
    for (const auto &cluster: all_clusters) {
        size_hist[cluster.Size()]++;
    }
    std::ofstream sizes(fs::append_path(output_dir, "cluster_sizes.tsv"));
    for (const auto &entry: size_hist) {
        sizes << entry.first << "\t" << entry.second << std::endl;
    }
    auto path_clusters = path_extractor_helper.GetPathClusters(all_clusters);
    INFO(path_clusters.size() << " path clusters");
    size_hist.clear();
    for (const auto &cluster: path_clusters) {
        size_hist[cluster.Size()]++;
    }
    std::ofstream path_sizes(fs::append_path(output_dir, "path_sizes.tsv"));
    for (const auto &entry: size_hist) {
        path_sizes << entry.first << "\t" << entry.second << std::endl;
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
                           const std::string &output_path) {
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

    std::ofstream os(fs::append_path(output_path, "normalized_tellseq_links.tsv"));
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

    fs::make_dir(cfg.output_dir);
    fs::make_dir(cfg.tmpdir);

    INFO("Loading graph");
    std::unique_ptr<io::IdMapper<std::string>> struct_id_mapper(new io::IdMapper<std::string>());
    std::unique_ptr<io::IdMapper<std::string>> id_mapper(new io::IdMapper<std::string>());
    //fixme optional multiplex reading
    cont_index::MultiplexGFAReader demulti_gfa(cfg.graph);
    gfa::GFAReader gfa(cfg.graph);
    INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links() << ", paths: " << gfa.num_paths());
    VERIFY_MSG(demulti_gfa.k() != -1U, "Failed to determine k-mer length");
    INFO(demulti_gfa.k());
//    VERIFY_MSG(gfa.k() % 2 == 1, "k-mer length must be odd");
    debruijn_graph::Graph struct_graph(demulti_gfa.k());
    gfa.to_graph(struct_graph, struct_id_mapper.get());
    debruijn_graph::Graph graph(demulti_gfa.k());
    demulti_gfa.to_graph(struct_graph, struct_id_mapper.get(), graph, id_mapper.get());
    INFO("Graph loaded. Total vertices: " << graph.size() << ", total edges: " << graph.e_size());

    std::ofstream graph_out(fs::append_path(cfg.output_dir, "internal_assembly_graph.gfa"));
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
        const size_t read_linkage_distance = 5000;

        //graph construction
        const double graph_score_threshold = 3.99;
        const size_t tail_threshold = 20000;
        const size_t length_threshold = 0;
        const size_t count_threshold = 3;

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

        std::unique_ptr<TimeTracerRAII> traceraii;
        traceraii.reset(new TimeTracerRAII(argv[0], 500));
        INFO("Time tracing is enabled");

        TIME_TRACE_SCOPE("Containment index");

//        GFAGraphConstructor gfa_graph_constructor(graph, gfa, id_mapper.get());
//        auto hifi_graph = gfa_graph_constructor.ConstructGraph();

        auto &lib = dataset[cfg.libindex];
        if (lib.type() == io::LibraryType::Clouds10x) {
            cont_index::ConstructBarcodeIndex(barcode_index,
                                              lib,
                                              graph,
                                              cfg.tmpdir,
                                              cfg.nthreads,
                                              frame_size,
                                              cfg.bin_load,
                                              cfg.debug);
        } else {
            WARN("Only read cloud libraries with barcode tags are supported for links");
        }
        INFO("Barcode index size: " << barcode_index.size());
        using BarcodeExtractor = barcode_index::FrameBarcodeIndexInfoExtractor;
        auto barcode_extractor_ptr = std::make_shared<BarcodeExtractor>(barcode_index, graph);

        size_t total_reads = 0;
        for (const auto &edge: graph.edges()) {
            auto begin = barcode_extractor_ptr->barcode_iterator_begin(edge);
            auto end = barcode_extractor_ptr->barcode_iterator_end(edge);
            for (auto it = begin; it != end; ++it) {
                total_reads += it->second.GetCount();
            }
        }
        INFO(total_reads << " total reads in barcode index");

//        scaffold_graph::ScaffoldGraph tellseq_graph(graph);

        //fixme move to separate function
//        std::unordered_map<std::string, EdgeId> seg_to_edge;
//        for (const EdgeId &edge: graph.canonical_edges()) {
//            seg_to_edge.emplace((*id_mapper)[edge.int_id()], edge);
//            seg_to_edge.emplace((*id_mapper)[edge.int_id()] + "'", graph.conjugate(edge));
//        }
//        std::unordered_set<std::pair<debruijn_graph::EdgeId, debruijn_graph::EdgeId>> correct_transitions;
//        std::vector<std::vector<debruijn_graph::EdgeId>> reference_paths;
//        std::ifstream ref_stream(cfg.refpath);
//        std::string path_string;
//        std::string ref_name;
//        size_t path_length;
//        std::string edge_name;
//        while (ref_stream >> path_string) {
//            std::vector<debruijn_graph::EdgeId> reference_path;
//            ref_name = path_string.substr(0, path_string.find(":"));
//            INFO(ref_name);
//            ref_stream >> path_length;
//            INFO(path_length);
//            for (size_t i = 0; i < path_length; ++i) {
//                ref_stream >> edge_name;
////                INFO(edge_name);
//                EdgeId edge = seg_to_edge.at(edge_name);
////                INFO(edge.int_id());
//                reference_path.push_back(edge);
//            }
//            INFO("Read reference path");
//            reference_paths.push_back(reference_path);
//        }
//        INFO("Getting gap paths");
//        size_t total_pairs = 0;
//        size_t connected_pairs = 0;
//        for (const auto &path: reference_paths) {
//            if (path.size() < 2)
//                continue;
//            for (auto first = path.begin(), second = std::next(path.begin()); second != path.end(); ++first, ++second) {
//                ++total_pairs;
//                EdgeId first_edge = *first;
//                EdgeId second_edge = *second;
////                INFO(first_edge.int_id() << ", " << second_edge.int_id());
//                if (graph.EdgeEnd(first_edge) == graph.EdgeEnd(second_edge)) {
//                    ++connected_pairs;
//                }
//            }
//        }
//        INFO(connected_pairs << " connected pairs out of " << total_pairs);
//move to separate function

        auto tellseq_graph = cont_index::GetTellSeqScaffoldGraph(graph, barcode_extractor_ptr, graph_score_threshold,
                                                                 length_threshold,
                                                                 tail_threshold, count_threshold, cfg.nthreads,
                                                                 cfg.bin_load,
                                                                 cfg.debug, cfg.output_dir, id_mapper.get());

//        LinkIndexGraphConstructor link_index_constructor(graph,
//                                                         barcode_extractor_ptr,
//                                                         graph_score_threshold,
//                                                         tail_threshold,
//                                                         length_threshold,
//                                                         count_threshold,
//                                                         cfg.nthreads);
//        auto score_function = link_index_constructor.ConstructScoreFunction();
        NormalizeTellseqLinks(tellseq_graph,
                              length_threshold,
                              barcode_extractor_ptr,
                              count_threshold,
                              tail_threshold,
                              id_mapper.get(),
                              cfg.output_dir);

//        GetLongEdgeStatistics(graph, barcode_index, training_length_threshold, training_length_offset,
//                              training_min_read_threshold, training_linkage_distance, cfg.nthreads, cfg.output_dir);
//        auto compare_output_path = fs::append_path(cfg.output_dir, "hifi_tellseq_scores.tsv");
//        CompareLinks(hifi_graph, tellseq_graph, score_function, id_mapper.get(), compare_output_path);

//        GetPathClusterStatistics(graph,
//                                 barcode_extractor_ptr,
//                                 hifi_graph,
//                                 read_linkage_distance,
//                                 relative_score_threshold,
//                                 min_read_threshold,
//                                 length_threshold,
//                                 cfg.nthreads, cfg.output_dir);


    }
}
