#pragma once

#include <common/io/reads/osequencestream.hpp>
#include <common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction_pipeline.hpp>
#include "common/assembly_graph/contracted_graph/contracted_graph_builder.hpp"
#include "statistics_processor.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph_constructor.hpp"
#include "modules/path_extend/read_cloud_path_extend/read_cloud_connection_conditions.hpp"

namespace scaffold_graph_utils {
    using path_extend::scaffold_graph::ScaffoldGraph;
    using path_extend::scaffold_graph::ScaffoldVertex;
    class ScaffoldGraphConstructor {
        const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
        size_t distance_;
        const Graph &g_;
        size_t max_threads_;

     public:
        ScaffoldGraphConstructor(const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage_, size_t distance_,
                                 const Graph &g_, size_t max_threads) : unique_storage_(unique_storage_), distance_(distance_), g_(g_),
                                                    max_threads_(max_threads) {}

        ScaffoldGraph ConstructScaffoldGraphUsingDijkstra() {
            ScaffoldGraph scaffold_graph(g_);
    //        auto bounded_dij = DijkstraHelper<Graph>::CreateBoundedDijkstra(g_, distance_, 10000);

            set<ScaffoldVertex> target_vertices;
            std::copy(unique_storage_.begin(), unique_storage_.end(), std::inserter(target_vertices, target_vertices.end()));

            path_extend::scaffold_graph::UniqueScaffoldGraphConstructor constructor(g_, unique_storage_,
                                                                                    target_vertices,
                                                                                    distance_, max_threads_);

//            for (const auto unique_edge: unique_storage_) {
//                scaffold_graph.AddVertex(unique_edge);
//            }
//            for (const auto unique_edge: unique_storage_) {
//                auto dij = omnigraph::CreateUniqueDijkstra(g_, distance_, unique_storage_);
//                dij.Run(g_.EdgeEnd(unique_edge));
//                for (auto v: dij.ReachedVertices()) {
//                    size_t distance = dij.GetDistance(v);
//                    if (distance < distance_) {
//                        for (auto connected: g_.OutgoingEdges(v)) {
//                            if (unique_storage_.IsUnique(connected) and connected != unique_edge
//                                and connected != g_.conjugate(unique_edge)) {
//                                ScaffoldGraph::ScaffoldEdge scaffold_edge(unique_edge, connected, (size_t) -1, 0, distance);
//                                scaffold_graph.AddEdge(scaffold_edge);
//                            }
//                        }
//                    }
//                }
//            }
            return *(constructor.Construct());
        }

        ScaffoldGraph ConstructScaffoldGraphFromContractedGraph(const contracted_graph::ContractedGraph &contracted_graph) {
            ScaffoldGraph scaffold_graph(g_);
            unordered_set<ScaffoldVertex> incoming_set;
            unordered_set<ScaffoldVertex> outcoming_set;
            unordered_set<ScaffoldVertex> union_set;
            unordered_set<VertexId> vertices;

            for (const auto &vertex: contracted_graph) {
                DEBUG("Vertex: " << vertex.int_id());
                vertices.insert(vertex);
                vector<ScaffoldVertex> incoming_vector;
                vector<ScaffoldVertex> outcoming_vector;
                for (auto it_in = contracted_graph.in_begin(vertex);
                     it_in != contracted_graph.in_end(vertex); ++it_in) {
                    for (const auto& edge: (*it_in).second) {
                        DEBUG("Incoming: " << edge.int_id());
                        incoming_vector.push_back(edge);
                        incoming_set.insert(edge);
                        union_set.insert(edge);
                        scaffold_graph.AddVertex(edge);
                    }
                }
                for (auto it_out = contracted_graph.out_begin(vertex);
                     it_out != contracted_graph.out_end(vertex); ++it_out) {
                    for (const auto& edge: (*it_out).second) {
                        DEBUG("Outcoming: " << edge.int_id());
                        outcoming_vector.push_back(edge);
                        outcoming_set.insert(edge);
                        union_set.insert(edge);
                        scaffold_graph.AddVertex(edge);
                    }
                }
                for (const auto& in_edge: incoming_vector) {
                    for (const auto& out_edge: outcoming_vector) {
                        DEBUG("Adding edge");
                        ScaffoldGraph::ScaffoldEdge scaffold_edge(in_edge, out_edge);
                        scaffold_graph.AddEdge(scaffold_edge);
                    }
                }
            }
            DEBUG("Incoming: " << incoming_set.size());
            DEBUG("Outcoming: " << outcoming_set.size());
            DEBUG("Union: " << union_set.size());
            return scaffold_graph;
        }

        ScaffoldGraph ConstructBarcodeScoreScaffoldGraph(const ScaffoldGraph& scaffold_graph,
                                                         shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> extractor,
                                                         const Graph& graph, size_t count_threshold,
                                                         size_t tail_threshold, double score_threshold) {
            auto score_function = make_shared<path_extend::NormalizedBarcodeScoreFunction>(graph, extractor,
                                                                                           count_threshold, tail_threshold);
            size_t num_threads = cfg::get().max_threads;
            INFO("Constructing barcode score scaffold graph, score threshold: " << score_threshold);
            path_extend::scaffold_graph::ScoreFunctionScaffoldGraphConstructor
                scaffold_graph_constructor(graph, scaffold_graph, score_function, score_threshold, num_threads);
            return *(scaffold_graph_constructor.Construct());
        }

        ScaffoldGraph ConstructLongGapScaffoldGraph(const ScaffoldGraph& scaffold_graph,
                                                    const path_extend::ScaffoldingUniqueEdgeStorage& storage,
                                                    shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> short_edge_extractor,
                                                    shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor,
                                                    const Graph& graph, const path_extend::ReadCloudMiddleDijkstraParams& params) {
            auto predicate = make_shared<path_extend::ReadCloudMiddleDijkstraPredicate>(graph, storage,
                                                                                        short_edge_extractor,
                                                                                        long_edge_extractor, params);
            size_t max_threads = cfg::get().max_threads;
            path_extend::scaffold_graph::PredicateScaffoldGraphConstructor constructor(graph, scaffold_graph,
                                                                                       predicate, max_threads);
            return *(constructor.Construct());
        }

        ScaffoldGraph ConstructCompositeConnectionScaffoldGraph(const ScaffoldGraph& scaffold_graph,
                                                                shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor,
                                                                shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor,
                                                                const path_extend::ScaffolderParams& params,
                                                                const debruijn_graph::conj_graph_pack& gp) {
            const size_t max_threads = cfg::get().max_threads;
            path_extend::CompositeConnectionConstructorCaller composite_constructor_caller(gp, main_extractor,
                                                                                           long_edge_extractor,
                                                                                           unique_storage_, max_threads);
            auto scaffold_graph_constructor = composite_constructor_caller.GetScaffoldGraphConstuctor(params, scaffold_graph);
            return *(scaffold_graph_constructor->Construct());
        }


        ScaffoldGraph ConstructOrderingScaffoldGraph(const ScaffoldGraph& scaffold_graph,
                                                     const barcode_index::FrameBarcodeIndexInfoExtractor& extractor,
                                                     const Graph& graph,
                                                     size_t count_threshold, double strictness) {
            auto predicate = make_shared<path_extend::EdgeSplitPredicate>(graph, extractor, count_threshold, strictness);
            size_t max_threads = cfg::get().max_threads;
            path_extend::scaffold_graph::PredicateScaffoldGraphConstructor constructor(graph, scaffold_graph,
                                                                                       predicate, max_threads);
            return *(constructor.Construct());
        }


        ScaffoldGraph ConstructNonTransitiveGraph(const ScaffoldGraph& scaffold_graph, const Graph& graph) {
            size_t distance_threshold = 3;
            size_t max_threads = 1;
            auto predicate = make_shared<path_extend::TransitiveEdgesPredicate>(scaffold_graph, graph, distance_threshold);
            path_extend::scaffold_graph::PredicateScaffoldGraphConstructor constructor(graph, scaffold_graph,
                                                                                       predicate, max_threads);
            return *(constructor.Construct());
        }

        DECL_LOGGER("ScaffoldGraphConstructor");
    };

    class OutDegreeDistribuiton: public read_cloud_statistics::Statistic {
        std::map<size_t, size_t> degree_distribution_;

     public:
        OutDegreeDistribuiton(): read_cloud_statistics::Statistic("out_degree_distribution"), degree_distribution_() {}
        OutDegreeDistribuiton(const OutDegreeDistribuiton& other) = default;
        void Insert(size_t degree) {
            degree_distribution_[degree]++;
        }

        void Serialize(const string& path) override {
            ofstream fout(path);
            for (const auto& entry: degree_distribution_) {
                fout << entry.first << " " << entry.second << std::endl;
            }
        }
    };

    class EdgeToDegree: public read_cloud_statistics::Statistic {
        typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;

        std::map<ScaffoldVertex, size_t> edge_to_degree_;

     public:
        EdgeToDegree(): read_cloud_statistics::Statistic("edge_to_degree"), edge_to_degree_() {}
        EdgeToDegree(const EdgeToDegree& other) = default;
        void Insert(const ScaffoldVertex& edge, size_t degree) {
            edge_to_degree_[edge] = degree;
        }

        void Serialize(const string& path) override {
            ofstream fout(path);
            for (const auto& entry: edge_to_degree_) {
                fout << entry.first.int_id() << " " << entry.second << std::endl;
            }
        }
    };

    class ScaffoldGraphAnalyzer : public read_cloud_statistics::StatisticProcessor {
        ScaffoldGraph graph_;

     public:
        ScaffoldGraphAnalyzer(const ScaffoldGraph& graph) : read_cloud_statistics::StatisticProcessor("scaffold_graph_stats"),
                                                            graph_(graph) {}

        void FillStatistics() override {
            auto out_distribution_ptr = std::make_shared<OutDegreeDistribuiton>(GetOutDegreeDistribution());
            auto edge_to_degree_ptr = std::make_shared<EdgeToDegree>(GetEdgeToDegree());
            AddStatistic(edge_to_degree_ptr);
            AddStatistic(out_distribution_ptr);
        }

        OutDegreeDistribuiton GetOutDegreeDistribution() const {
            OutDegreeDistribuiton distribution;
            for (const ScaffoldGraph::ScaffoldGraphVertex &vertex: graph_.vertices()) {
                distribution.Insert(graph_.OutgoingEdgeCount(vertex));
            }
            return distribution;
        }

        EdgeToDegree GetEdgeToDegree() const {
            EdgeToDegree edge_to_degree;
            for (const ScaffoldGraph::ScaffoldGraphVertex& vertex: graph_.vertices()) {
                size_t degree = graph_.OutgoingEdgeCount(vertex);
                edge_to_degree.Insert(vertex, degree);
            }
            return edge_to_degree;
        }

        DECL_LOGGER("ScaffoldGraphAnalyzer")
    };

    struct SignedScaffoldGraph {
      typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
      ScaffoldGraph graph_;
      std::unordered_map<ScaffoldVertex, bool> vertex_to_sign_;

      SignedScaffoldGraph(const ScaffoldGraph& graph_, const unordered_map<ScaffoldVertex, bool>& vertex_to_sign_) : graph_(
          graph_), vertex_to_sign_(vertex_to_sign_) {}
    };

    class ScaffoldGraphSplitter {
     public:
        typedef path_extend::validation::EdgeWithMapping EdgeWithMapping;
        typedef vector<EdgeWithMapping> ContigPath;
        typedef path_extend::validation::NamedSimplePath NamedSimplePath;
        typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
     private:
        const Graph& g_;
        const vector<NamedSimplePath>& reference_paths_;

     public:
        ScaffoldGraphSplitter(const Graph& g_, const vector<NamedSimplePath>& named_reference_paths)
            : g_(g_), reference_paths_(named_reference_paths) {}

        map<string, SignedScaffoldGraph> SplitAndSignScaffoldGraph(const ScaffoldGraph& graph) {
            map<string, SignedScaffoldGraph> result;
            for (const auto& named_path: reference_paths_) {
                string name = named_path.name_;
                ContigPath path = named_path.path_;
                const size_t MIN_SIZE = 10;
                INFO("Name: " << name);
                INFO("Path size: " << path.size());
                if (not ContainsPathOrRC(result, name) and path.size() >= MIN_SIZE) {
                    auto path_edges = ExtractEdges(path);
                    ScaffoldGraph subgraph(g_);
                    unordered_map<ScaffoldVertex, bool> vertex_to_sign;
                    for (const auto& vertex: graph.vertices()) {
                        if (IsPlusVertex(path_edges, vertex)) {
                            vertex_to_sign.insert({vertex, true});
                        }
                        if (IsMinusVertex(path_edges, vertex)) {
                            VERIFY(vertex_to_sign.find(vertex) == vertex_to_sign.end());
                            vertex_to_sign.insert({vertex, false});
                        }
                        if (vertex_to_sign.find(vertex) != vertex_to_sign.end()) {
                            subgraph.AddVertex(vertex);
                        }
                    }
                    for (const ScaffoldGraph::ScaffoldEdge& edge: graph.edges()) {
                        if (subgraph.Exists(edge.getStart()) and subgraph.Exists(edge.getEnd())) {
                            subgraph.AddEdge(edge);
                        }
                    }
                    SignedScaffoldGraph signed_subgraph(subgraph, vertex_to_sign);
                    result.insert({name, signed_subgraph});
                }
            }
            return result;
        }

        bool ContainsPathOrRC(const map<string, SignedScaffoldGraph>& name_to_graph, const string& name) {
            const string rc_suffix = "_RC";
            string stripped_name = name.substr(0, name.size() - rc_suffix.size());

            return name_to_graph.find(name + rc_suffix) != name_to_graph.end() or
                name_to_graph.find(stripped_name) != name_to_graph.end();
        }

        unordered_set<ScaffoldVertex> ExtractEdges(const ContigPath& simple_path) {
            unordered_set<ScaffoldVertex> result;
            for (const auto& ewm: simple_path) {
                result.insert(ewm.edge_);
            }
            return result;
        }

        bool IsPlusVertex(const std::unordered_set<ScaffoldVertex>& path_edges, const ScaffoldVertex& vertex) {
            return path_edges.find(vertex) != path_edges.end();
        }

        bool IsMinusVertex(const std::unordered_set<ScaffoldVertex>& path_edges, const ScaffoldVertex& vertex) {
            return path_edges.find(vertex.getConjugateFromGraph(g_)) != path_edges.end();
        }
    };

    class ScaffoldGraphPrinter {
     public:
        typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
        typedef path_extend::validation::EdgeWithMapping EdgeWithMapping;
        typedef vector<EdgeWithMapping> ContigPath;
     private:
        const Graph& g_;
     public:
        ScaffoldGraphPrinter(const Graph& g_)
            : g_(g_) {}
     public:
        void PrintGraph(const SignedScaffoldGraph& signed_graph, const string& contig_path, const string& graph_path) {
            auto vertex_to_name = GetVertexToName(signed_graph.graph_);
            PrintGraphAsContigs(signed_graph.graph_, contig_path, vertex_to_name);
            PrintGraphAsGraph(signed_graph, graph_path, vertex_to_name);
        }

        void PrintMultipleGraphs(const map<string, SignedScaffoldGraph>& name_to_signed_graph, const string& output_path) {
            string contig_path = fs::append_path(output_path, "contigs");
            string graph_path = fs::append_path(output_path, "graphs");
            fs::make_dir(contig_path);
            fs::make_dir(graph_path);
            for (const auto& entry: name_to_signed_graph) {
                string contig_name = fs::append_path(contig_path, entry.first);
                string graph_name = fs::append_path(graph_path, entry.first);
                PrintGraph(entry.second, contig_name, graph_name);
            }
        }

     private:
        void PrintGraphAsContigs(const ScaffoldGraph& graph, const string& filename,
                                 const std::unordered_map<ScaffoldVertex, string>& vertex_to_name) {
            ofstream fout(filename + ".fasta");

            for (const ScaffoldVertex& vertex: graph.vertices()) {
                fout << ">" << vertex_to_name.at(vertex) << endl;
                path_extend::scaffold_graph::EdgeGetter getter;
                EdgeId edge = getter.GetEdgeFromScaffoldVertex(vertex);
                io::WriteWrapped(g_.EdgeNucls(edge).str(), fout);
            }
        }


        void PrintGraphAsGraph(const SignedScaffoldGraph& signed_graph, const string& filename,
                               const std::unordered_map<ScaffoldVertex, string>& vertex_to_name) {
            ofstream fout(filename + ".scg");
            for (const ScaffoldGraph::ScaffoldEdge& edge: signed_graph.graph_.edges()) {
                ScaffoldVertex start = edge.getStart();
                ScaffoldVertex end = edge.getEnd();
                string start_sign = signed_graph.vertex_to_sign_.at(start) ? "+" : "-";
                string end_sign = signed_graph.vertex_to_sign_.at(end) ? "+" : "-";
                fout << "(" << vertex_to_name.at(start) << "_" << start_sign << ") ";
                fout << "(" << vertex_to_name.at(end) << "_" << end_sign << ") ";
                fout << edge.getWeight() << " " << edge.getLength() << endl;
            }
        }

        std::unordered_map<ScaffoldVertex, string> GetVertexToName(const ScaffoldGraph& graph) {
            std::unordered_map<ScaffoldVertex, string> vertex_to_name;
            size_t id = 0;
            for (const ScaffoldVertex& vertex: graph.vertices()) {
                string name = io::MakeContigId(id++, vertex.getLengthFromGraph(g_), vertex.getCoverageFromGraph(g_));
                vertex_to_name.insert({vertex, name});
            }
            return vertex_to_name;
        };
    };

    class MetaScaffoldGraphPrinter {
        const Graph& g_;

     public:
        MetaScaffoldGraphPrinter(const Graph& g_) : g_(g_) {}

        void PrintGraphAsMultiple(const ScaffoldGraph& graph, const vector<path_extend::validation::NamedSimplePath>& reference_paths,
                             const string& output_path) {
            ScaffoldGraphSplitter splitter(g_, reference_paths);
            auto name_to_signed_graph = splitter.SplitAndSignScaffoldGraph(graph);
            ScaffoldGraphPrinter printer(g_);
//            printer.PrintMultipleGraphs(name_to_signed_graph, output_path);
        }
    };
}