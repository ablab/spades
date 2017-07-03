#pragma once
#include <boost/pending/disjoint_sets.hpp>
#include "assembly_graph/graph_support/scaff_supplementary.hpp"
#include "transitions.hpp"
#include "statistics_processor.hpp"

namespace contracted_graph {

    class AdjacencyMap {
        std::map<VertexId, vector<EdgeId>> data_;

    public:
        typedef std::map<VertexId, vector<EdgeId>>::const_iterator const_iterator;
        AdjacencyMap() = default;
        AdjacencyMap(const VertexId& vertex, const EdgeId& edge) : data_({{vertex, {edge}}}) {}
        void InsertPair(const VertexId& vertex, const EdgeId& edge) {
            data_[vertex].push_back(edge);
        }

        const_iterator begin() const {
            return data_.begin();
        }

        const_iterator end() const {
            return data_.end();
        }

    };

    //todo: replace with boost graph?
    class ContractedGraph {
        std::map<VertexId, AdjacencyMap> outcoming_;
        std::map<VertexId, AdjacencyMap> incoming_;
    public:

        typedef std::map<VertexId, AdjacencyMap>::const_iterator const_iterator;

        ContractedGraph() = default;
        void InsertEdge(const VertexId& head, const VertexId& tail, const EdgeId& edge) {
            outcoming_[head].InsertPair(tail, edge);
            incoming_[tail].InsertPair(head, edge);
            AdjacencyMap empty;
            if (outcoming_.find(tail) == outcoming_.end()) {
                outcoming_[tail] = empty;
            }
            if (incoming_.find(head) == incoming_.end()) {
                incoming_[head] = empty;
            }
        }

        AdjacencyMap::const_iterator incoming_begin(const VertexId& vertex) const {
            return incoming_.at(vertex).begin();
        }

        AdjacencyMap::const_iterator incoming_end(const VertexId& vertex) const {
            return incoming_.at(vertex).end();
        }

        AdjacencyMap::const_iterator outcoming_begin(const VertexId& vertex) const {
            return outcoming_.at(vertex).begin();
        }

        AdjacencyMap::const_iterator outcoming_end(const VertexId& vertex) const {
            return outcoming_.at(vertex).end();
        }

        vector <EdgeId> GetIncoming(const VertexId& vertex) {
            vector<EdgeId> incoming;
            for (auto in_it = incoming_begin(vertex); in_it != incoming_end(vertex); ++in_it) {
                for (auto edge_it = (*in_it).second.begin(); edge_it != (*in_it).second.end(); ++edge_it) {
                    incoming.push_back(*edge_it);
                }
            }
            return incoming;
        }

        vector <EdgeId> GetOutcoming(const VertexId& vertex) {
            vector<EdgeId> outcoming;
            for (auto out_it = outcoming_begin(vertex); out_it != outcoming_end(vertex); ++out_it) {
                for (auto edge_it = (*out_it).second.begin(); edge_it != (*out_it).second.end(); ++edge_it) {
                    outcoming.push_back(*edge_it);
                }
            }
            return outcoming;
        }

        const_iterator begin() const {
            return outcoming_.begin();
        }

        const_iterator end() const {
            return outcoming_.end();
        }
    };

    class ContractedGraphBuilder {
        const Graph &g_;
        path_extend::ScaffoldingUniqueEdgeStorage unique_storage_;
        vector<EdgeId> long_edges_;
        unordered_set<VertexId> long_vertices_;

    public:
        typedef std::map<VertexId, std::size_t> rank_t; // => order on Element
        typedef std::map<VertexId, VertexId> parent_t;
        typedef boost::associative_property_map<rank_t> rank_property_t;
        typedef boost::associative_property_map<parent_t> parent_property_t;
        typedef boost::disjoint_sets<rank_property_t, parent_property_t> dsu_t;


        ContractedGraphBuilder(const Graph &g,
                               const path_extend::ScaffoldingUniqueEdgeStorage unique_storage) : g_(g),
                                                                                                 unique_storage_(unique_storage) {}

        ContractedGraph BuildContractedGraph() {
            rank_t rank_map;
            parent_t parent_map;
            rank_property_t rank_pmap(rank_map);
            parent_property_t parent_pmap(parent_map);
            omnigraph::IterationHelper<Graph, VertexId> vertex_iteration_helper(g_);

            auto dsu = BuildDSU(g_, rank_pmap, parent_pmap);

            ContractedGraph graph;

            GetConnectedGraphStats(dsu);
            GetLongEdgesStats(long_edges_);

            for (const auto& edge: long_edges_) {
                graph.InsertEdge(dsu.find_set(g_.EdgeStart(edge)), dsu.find_set(g_.EdgeEnd(edge)), edge);
            }

            return graph;
        }

    private:

        void GetLongEdgesStats(const vector<EdgeId>& long_edges) {
            size_t total_length = 0;
            vector<size_t> length_distribution;
            for (const auto& edge: long_edges) {
                total_length += g_.length(edge);
                length_distribution.push_back(g_.length(edge));
            }
            std::sort(length_distribution.begin(), length_distribution.end());
            size_t n50 = GetN50(length_distribution, total_length);
            INFO("Long edge total length: " << total_length);
            INFO("Long edge N50: " << n50);
        }

        size_t GetN50(const vector<size_t>& length_distribution, const size_t total_length) {
            size_t current_length = 0;
            size_t half_length = total_length / 2;
            size_t result;
            for (auto it = length_distribution.begin(); current_length < half_length; ++it) {
                result = *it;
                current_length += *it;
            }
            return result;
        }

        void GetConnectedGraphStats(dsu_t &dsu) {
            omnigraph::IterationHelper<Graph, VertexId> vertex_iteration_helper(g_);
            omnigraph::IterationHelper<Graph, EdgeId> edge_iteration_helper(g_);
            unordered_set<VertexId> vertices;
            unordered_set<EdgeId> edges;

            for (auto it = vertex_iteration_helper.begin(); it != vertex_iteration_helper.end(); ++it) {
                VertexId vertex = *it;
                if (long_vertices_.find(dsu.find_set(vertex)) != long_vertices_.end()) {
                    vertices.insert(vertex);
                }
            }
            for (auto it = edge_iteration_helper.begin(); it != edge_iteration_helper.end(); ++it) {
                EdgeId edge = *it;
                VertexId start = dsu.find_set(g_.EdgeStart(edge));
                VertexId end = dsu.find_set(g_.EdgeEnd(edge));
                if (long_vertices_.find(start) != long_vertices_.end() or long_vertices_.find(end) != long_vertices_.end()) {
                    edges.insert(edge);
                }
            }
            vector<size_t> length_distribution;
            size_t total_length = 0;
            for (const auto& edge: edges) {
                length_distribution.push_back(g_.length(edge));
                total_length += g_.length(edge);
            }
            std::sort(length_distribution.begin(), length_distribution.end());

            size_t n50 = GetN50(length_distribution, total_length);
            INFO("True vertices: " << vertices.size());
            INFO("True edges: " << edges.size());
            INFO("True total length: " << total_length);
            INFO("N50: " << n50);
        }

        vector<EdgeId> GetLongEdges() const {
            return long_edges_;
        }

        template<typename Rank, typename Parent>
        boost::disjoint_sets<Rank, Parent> BuildDSU(const Graph &g, Rank &r, Parent &p) {
            omnigraph::IterationHelper<Graph, EdgeId> edge_iteration_helper(g);

            boost::disjoint_sets<Rank, Parent> dsets(r, p);
            for (auto it = g.vertices().begin(); it != g.vertices().end(); ++it) {
                dsets.make_set(*it);
            }
            INFO("Overall vertices: " << dsets.count_sets(g.vertices().begin(), g.vertices().end()));
            size_t short_edges = 0;
            size_t long_edges = 0;
            size_t self_linkages = 0;
            size_t loops = 0;
            unordered_set<VertexId> long_roots;
            for (auto it = edge_iteration_helper.begin(); it != edge_iteration_helper.end(); ++it) {
                auto edge = *it;
                if (not unique_storage_.IsUnique(edge)) {
                    ++short_edges;
                    VertexId start = g.EdgeStart(edge);
                    VertexId end = g.EdgeEnd(edge);
                    if (dsets.find_set(start) == dsets.find_set(end)) {
                        ++self_linkages;
                    }
                    dsets.union_set(g.EdgeStart(edge), g.EdgeEnd(edge));
                } else {
                    ++long_edges;
                    if (dsets.find_set(g.EdgeStart(edge)) == dsets.find_set(g.EdgeEnd(edge))) {
                        loops++;
                    }
                    long_roots.insert(dsets.find_set(g.EdgeEnd(edge)));
                    long_roots.insert(dsets.find_set(g.EdgeStart(edge)));
                    long_edges_.push_back(edge);
                }
            }
            long_vertices_ = long_roots;
            INFO("CONTRACTED GRAPH STATISTICS");
//            INFO("Short edges: " << short_edges);
//            INFO("Self linkages: " << self_linkages);
//            INFO("Vertices in contracted: " << dsets.count_sets(g.vertices().begin(), g.vertices().end()));
            INFO("Vertices: " << long_roots.size());
            INFO("Edges: " << long_edges);
            INFO("Loops: " << loops);
            return dsets;
        };
    };
class OutDegreeDistribuiton: public read_cloud_statistics::Statistic {
    std::map<size_t, size_t> degree_distribution_;

 public:
    OutDegreeDistribuiton(): read_cloud_statistics::Statistic("out_degree_distribution"), degree_distribution_() {}
    OutDegreeDistribuiton(const OutDegreeDistribuiton& other) = default;
    void Insert(size_t degree) {
        degree_distribution_[degree]++;
    }

    void Serialize(ofstream& fout) override {
        for (const auto& entry: degree_distribution_) {
            fout << entry.first << " " << entry.second << std::endl;
        }
    }
};

}