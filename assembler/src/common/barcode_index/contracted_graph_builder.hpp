#pragma once

#include "contracted_graph.hpp"
#include "adt/concurrent_dsu.hpp"
#include <common/barcode_index/cluster_storage.hpp>
#include "assembly_graph/graph_support/scaff_supplementary.hpp"

namespace contracted_graph {

class ContractedGraphFactory {
 protected:
    shared_ptr<ContractedGraph> graph_ptr_;
 public:
    ContractedGraphFactory() : graph_ptr_(make_shared<ContractedGraph>()) {}
    virtual ~ContractedGraphFactory() = default;
    virtual void Construct() = 0;
    shared_ptr<ContractedGraph> GetGraph() {
        return graph_ptr_;
    }
};

class PartsBasedContractedFactory : public ContractedGraphFactory{

 protected:
    const Graph& g_;
    using ContractedGraphFactory::graph_ptr_;

    struct ContractedGraphParts {
      vector <EdgeId> long_edges_;
      unordered_set <VertexId> long_edge_ends_;
      unordered_map <VertexId, size_t> vertex_to_capacity_;
      unordered_map <VertexId, VertexId> vertex_to_root_;
    };

    virtual ContractedGraphParts ConstructParts() const = 0;

    void ConstructFromParts(ContractedGraphParts&& parts) {
        DEBUG("Constructing from parts");
        const auto& vertex_to_root = parts.vertex_to_root_;
        const auto& long_edges = parts.long_edges_;
        const auto& vertex_to_capacity = parts.vertex_to_capacity_;
        DEBUG("Vertex to root size: " << vertex_to_root.size());
        DEBUG("Vertex to capacity: " << vertex_to_capacity.size())
        for (const auto& edge: long_edges) {
            DEBUG("Processing edge " << edge.int_id());
            DEBUG(g_.EdgeStart(edge) << " -> " << g_.EdgeEnd(edge));
            VertexId start_root = vertex_to_root.at(g_.EdgeStart(edge));
            VertexId end_root = vertex_to_root.at(g_.EdgeEnd(edge));
            DEBUG("Inserting vertices and edges");
            this->graph_ptr_->InsertVertex(start_root);
            this->graph_ptr_->InsertVertex(end_root);
            this->graph_ptr_->InsertEdge(start_root, end_root, edge);
            this->graph_ptr_->InsertCapacity(start_root, vertex_to_capacity.at(start_root));
            this->graph_ptr_->InsertCapacity(end_root, vertex_to_capacity.at(end_root));
        }
    }
    DECL_LOGGER("DSUBasedContractedGraphFactory");
 public:
    PartsBasedContractedFactory(const Graph& assembly_graph_) : g_(assembly_graph_) {}
    virtual ~PartsBasedContractedFactory() {}

    void Construct() override {
        ConstructFromParts(ConstructParts());
    }
};

class DBGContractedGraphFactory : public PartsBasedContractedFactory {
    using PartsBasedContractedFactory::g_;
    using PartsBasedContractedFactory::graph_ptr_;
    using PartsBasedContractedFactory::ContractedGraphParts;
    typedef dsu::ConcurrentDSU contracted_dsu_t;

    const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;

 public:

    DBGContractedGraphFactory(const Graph& g,
        const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage) : PartsBasedContractedFactory(g), unique_storage_(unique_storage) {}

 private:
    ContractedGraphParts ConstructParts() const override {
        omnigraph::IterationHelper<Graph, EdgeId> edge_iteration_helper(g_);
        omnigraph::IterationHelper<Graph, VertexId> vertex_iteration_helper(g_);

        DEBUG("Preparing parts");
        ContractedGraphParts graph_parts;
        auto& unique_edge_ends = graph_parts.long_edge_ends_;
        auto& vertex_to_capacity = graph_parts.vertex_to_capacity_;
        auto& unique_edges = graph_parts.long_edges_;
        auto& vertex_to_root = graph_parts.vertex_to_root_;

        contracted_dsu_t graph_dsu(g_.size());
        std::unordered_map<VertexId, size_t> vertex_to_id;
        std::unordered_map<size_t, VertexId> id_to_vertex;

        size_t counter = 0;
        for (auto vertex : vertex_iteration_helper) {
            vertex_to_id.insert({vertex, counter});
            id_to_vertex.insert({counter, vertex});
            vertex_to_capacity.insert({vertex, 0});
        }

        DEBUG("Filling parts");
        for (auto it = edge_iteration_helper.begin(); it != edge_iteration_helper.end(); ++it) {
            auto edge = *it;
            VertexId start = g_.EdgeStart(edge);
            VertexId end = g_.EdgeEnd(edge);
            size_t start_id = vertex_to_id.at(start);
            size_t end_id = vertex_to_id.at(end);
            size_t start_root = graph_dsu.find_set(start_id);
            size_t end_root = graph_dsu.find_set(end_id);
            VertexId start_root_vertex = id_to_vertex.at(start_root);
            VertexId end_root_vertex = id_to_vertex.at(start_root);

            if (not unique_storage_.IsUnique(edge)) {
                if (start_root == end_root) {
                    vertex_to_capacity.at(start_root_vertex) += g_.length(edge);
                } else {
                    size_t start_capacity = vertex_to_capacity[start_root_vertex];
                    size_t end_capacity = vertex_to_capacity[end_root_vertex];
                    graph_dsu.unite(start_root, end_root);
                    VertexId new_vertex = id_to_vertex.at(graph_dsu.find_set(start_root));
                    vertex_to_capacity.at(start_root_vertex) = start_capacity + end_capacity + g_.length(edge);
                }
            } else {
                unique_edge_ends.insert(start_root_vertex);
                unique_edge_ends.insert(end_root_vertex);
                unique_edges.push_back(edge);
            }
        }
        DEBUG(graph_dsu.num_sets() << " sets in dsu");
        for (auto vertex: vertex_iteration_helper) {
            VertexId root_vertex = id_to_vertex.at(graph_dsu.find_set(vertex_to_id.at(vertex)));
            DEBUG("Inserting vertex and root: " << vertex.int_id() << ", " << root_vertex.int_id());
            vertex_to_root.insert({vertex, root_vertex});
        }
        return graph_parts;
    }

    DECL_LOGGER("DBGContractedGraphFactory");
};

class TransposedContractedGraphFactory: public ContractedGraphFactory {
    const ContractedGraph& other_;
 public:
    explicit TransposedContractedGraphFactory(const ContractedGraph& other) : other_(other) {}

    void Construct() override {
        TransposeContractedGraph(other_);
    }

 private:
    void TransposeContractedGraph(const ContractedGraph& other) {
        for (const auto& vertex: other) {
            graph_ptr_->InsertVertex(vertex);
        }

        for (const auto& start: other) {
            for (auto it = other.out_begin(start); it != other.out_end(start); ++it) {
                VertexId end = (*it).first;
                for (const auto& edge: (*it).second) {
                    graph_ptr_->InsertEdge(end, start, edge);
                }
            }
        }
    }
};

class SubgraphContractedGraphFactory: public ContractedGraphFactory {
    const ContractedGraph& other_;
    const std::unordered_set<VertexId>& vertices_;
 public:
    SubgraphContractedGraphFactory(const ContractedGraph& other, const std::unordered_set<VertexId>& vertices) :
        other_(other), vertices_(vertices) {}

    void Construct() override {
        ExtractSubgraphFromContractedGraph(other_, vertices_);
    }

 private:
    void ExtractSubgraphFromContractedGraph(const ContractedGraph& other, const std::unordered_set<VertexId>& vertices) {
        for (const auto& vertex: vertices) {
            VERIFY(other.ContainsVertex(vertex));
            graph_ptr_->InsertVertex(vertex);
        }

        for (const auto& vertex: vertices) {
            for (const auto& adj_list: other.outcoming(vertex)) {
                VertexId next = adj_list.first;
                if (vertices.find(next) != vertices.end()) {
                    for (const auto& edge: adj_list.second) {
                        graph_ptr_->InsertEdge(vertex, next, edge);
                    }
                }
            }
        }
    };
};


} //contracted_graph