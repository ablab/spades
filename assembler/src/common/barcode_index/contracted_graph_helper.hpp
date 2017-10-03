#pragma once

#include "contracted_graph_builder.hpp"

namespace contracted_graph {
class SimpleContractedGraphFactory: public PartsBasedContractedFactory {
    using PartsBasedContractedFactory::g_;
    using PartsBasedContractedFactory::graph_ptr_;
    using PartsBasedContractedFactory::ContractedGraphParts;
    typedef dsu::ConcurrentDSU contracted_dsu_t;
    typedef cluster_storage::Cluster::InternalGraph InternalGraph;
    const InternalGraph& internal_graph_;

 public:
    SimpleContractedGraphFactory(const Graph& assembly_graph_, const InternalGraph& internal_graph_)
        : PartsBasedContractedFactory(assembly_graph_), internal_graph_(internal_graph_) {}

 private:
    ContractedGraphParts ConstructParts() const override {
        TRACE("Building contracted graph")
        ContractedGraphParts graph_parts;
        std::unordered_set <VertexId> vertices;
        auto& vertex_to_root = graph_parts.vertex_to_root_;
        auto& vertex_to_capacity = graph_parts.vertex_to_capacity_;
        auto& long_edges = graph_parts.long_edges_;

        TRACE("Internal graph:")
        for (const auto& start: internal_graph_) {
            for (auto it = internal_graph_.outcoming_begin(start); it != internal_graph_.outcoming_end(start); ++it) {
                auto end = *it;
                TRACE(start.int_id() << " -> " << end.int_id());
            }
        }
        contracted_dsu_t graph_dsu(internal_graph_.size());
        std::unordered_map<VertexId, size_t> vertex_to_id;
        std::unordered_map<size_t, VertexId> id_to_vertex;

        size_t counter = 0;


        for (auto it = internal_graph_.begin(); it != internal_graph_.end(); ++it) {
            EdgeId edge = *it;
            vertices.insert(g_.EdgeStart(edge));
            vertices.insert(g_.EdgeEnd(edge));
            vertex_to_capacity[g_.EdgeStart(edge)] = 0;
            vertex_to_capacity[g_.EdgeEnd(edge)] = 0;
            vertex_to_id.insert({g_.EdgeStart(edge), counter});
            id_to_vertex.insert({counter, g_.EdgeEnd(edge)});
            long_edges.push_back(edge);
        }

        unordered_set <VertexId> long_roots;
        for (const auto& start: internal_graph_) {
            for (auto it = internal_graph_.outcoming_begin(start); it != internal_graph_.outcoming_end(start); ++it) {
                auto end = *it;
                VertexId start_vertex = g_.EdgeEnd(start);
                VertexId end_vertex = g_.EdgeStart(end);
                size_t start_id = vertex_to_id.at(start_vertex);
                size_t end_id = vertex_to_id.at(end_vertex);
                size_t start_root = graph_dsu.find_set(start_id);
                size_t end_root = graph_dsu.find_set(end_id);

                graph_dsu.unite(start_root, end_root);
            }
        }

        for (const auto& vertex: vertices) {
            VertexId root = id_to_vertex.at(graph_dsu.find_set(vertex_to_id.at(vertex)));
            vertex_to_root.insert({vertex, root});
        }
        return graph_parts;
    }
    DECL_LOGGER("SimpleContractedGraphFactory");
};

class ContractedGraphFactoryHelper {
 private:
    const Graph& g_;
 public:
    explicit ContractedGraphFactoryHelper(const Graph& g_) : g_(g_) {}

    ContractedGraph ConstructFromUniqueStorage(const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage) const {
        DBGContractedGraphFactory factory(g_, unique_storage);
        factory.Construct();
        return *(factory.GetGraph());
    }

    ContractedGraph ConstructFromInternalGraph(const cluster_storage::Cluster::SimpleGraph<EdgeId>& internal_graph) const {
        SimpleContractedGraphFactory factory(g_, internal_graph);
        factory.Construct();
        return *(factory.GetGraph());
    }

    ContractedGraph TransposeContractedGraph(const ContractedGraph& other) const {
        TransposedContractedGraphFactory factory(other);
        factory.Construct();
        return *(factory.GetGraph());
    }

    ContractedGraph ExtractContractedSubgraph(const ContractedGraph& other,
                                              const std::unordered_set<VertexId>& vertices) const {
        SubgraphContractedGraphFactory factory(other, vertices);
        factory.Construct();
        return *(factory.GetGraph());
    }
};
}