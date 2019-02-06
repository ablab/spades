#include "contracted_graph_builder.hpp"

namespace contracted_graph {

void PartsBasedContractedFactory::ConstructFromParts(PartsBasedContractedFactory::ContractedGraphParts&& parts) {
    DEBUG("Constructing from parts");
    const auto& vertex_to_root = parts.vertex_to_root_;
    const auto& long_edges = parts.long_edges_;
    const auto& vertex_to_capacity = parts.vertex_to_capacity_;
    DEBUG("Vertex to root size: " << vertex_to_root.size());
    DEBUG("Vertex to capacity: " << vertex_to_capacity.size());
    for (const auto& entry: vertex_to_root) {
        this->graph_ptr_->InsertVertex(entry.second);
    }
    for (const auto& edge: long_edges) {
        DEBUG("Processing edge " << edge.int_id());
        DEBUG(edge.GetStartGraphVertex(g_) << " -> " << edge.GetEndGraphVertex(g_));
        VertexId start_root = vertex_to_root.at(edge.GetStartGraphVertex(g_));
        VertexId end_root = vertex_to_root.at(edge.GetEndGraphVertex(g_));
        DEBUG("Inserting vertices and edges");
        this->graph_ptr_->InsertEdge(start_root, end_root, edge);

    }
    for (const auto& entry: parts.vertex_to_capacity_) {
        this->graph_ptr_->InsertCapacity(entry.first, entry.second);
    }
}
void PartsBasedContractedFactory::Construct() {
    ConstructFromParts(ConstructParts());
}

PartsBasedContractedFactory::ContractedGraphParts SimpleContractedGraphFactory::ConstructParts() const {
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
    std::unordered_map<VertexId, size_t> vertex_to_id;
    std::unordered_map<size_t, VertexId> id_to_vertex;

    size_t counter = 0;

    for (auto it = internal_graph_.begin(); it != internal_graph_.end(); ++it) {
        SimpleContractedGraphFactory::ScaffoldVertex edge = *it;
        VertexId start = edge.GetStartGraphVertex(g_);
        VertexId end = edge.GetEndGraphVertex(g_);
        vertices.insert(start);
        vertices.insert(end);
        vertex_to_capacity[start] = 0;
        vertex_to_capacity[end] = 0;
        if (vertex_to_id.find(start) == vertex_to_id.end()) {
            vertex_to_id.insert({start, counter});
            id_to_vertex.insert({counter, start});
            ++counter;
            TRACE("Inserting vertex " << start.int_id() << " with index " << counter);
        }
        if (vertex_to_id.find(end) == vertex_to_id.end()) {
            vertex_to_id.insert({end, counter});
            id_to_vertex.insert({counter, end});
            ++counter;
            TRACE("Inserting vertex " << end.int_id() << " with index " << counter);
        }
        long_edges.push_back(edge);
    }

    VERIFY(counter > 0);
    contracted_dsu_t graph_dsu(counter);
    TRACE("Found " << counter << " vertices");
    DEBUG("Extracted vertices");

    unordered_set <VertexId> long_roots;
    for (const auto& start: internal_graph_) {
        for (auto it = internal_graph_.outcoming_begin(start); it != internal_graph_.outcoming_end(start); ++it) {
            auto end = *it;
            VertexId start_vertex = start.GetEndGraphVertex(g_);
            VertexId end_vertex = end.GetStartGraphVertex(g_);
            TRACE("Merging vertices " << start_vertex.int_id() << " and " << end_vertex.int_id());
            size_t start_id = vertex_to_id.at(start_vertex);
            size_t end_id = vertex_to_id.at(end_vertex);
            size_t start_root = graph_dsu.find_set(start_id);
            size_t end_root = graph_dsu.find_set(end_id);
            graph_dsu.unite(start_root, end_root);
        }
    }

    DEBUG("Built dsu");

    for (const auto& vertex: vertices) {
        VertexId root = id_to_vertex.at(graph_dsu.find_set(vertex_to_id.at(vertex)));
        vertex_to_root.insert({vertex, root});
    }
    return graph_parts;
}

PartsBasedContractedFactory::ContractedGraphParts DBGContractedGraphFactory::ConstructParts() const {
        omnigraph::IterationHelper<Graph, EdgeId> edge_iteration_helper(g_);
        omnigraph::IterationHelper<Graph, VertexId> vertex_iteration_helper(g_);
    DEBUG("Preparing parts");
    ContractedGraphParts graph_parts;

    contracted_dsu_t graph_dsu(g_.size());
    std::unordered_map<VertexId, size_t> vertex_to_id;
    std::unordered_map<size_t, VertexId> id_to_vertex;

    size_t counter = 0;
    for (auto vertex : vertex_iteration_helper) {
        vertex_to_id.insert({vertex, counter});
        id_to_vertex.insert({counter, vertex});
        graph_parts.vertex_to_capacity_.insert({vertex, 0});
        ++counter;
    }

    DEBUG("Filling parts");
    for (auto it = edge_iteration_helper.begin(); it != edge_iteration_helper.end(); ++it) {
        auto edge = *it;

        ProcessEdge(graph_dsu, graph_parts, vertex_to_id, id_to_vertex, edge);
    }
    DEBUG(graph_dsu.num_sets() << " sets in dsu");
    for (auto vertex: vertex_iteration_helper) {
        VertexId root_vertex = id_to_vertex.at(graph_dsu.find_set(vertex_to_id.at(vertex)));
        DEBUG("Inserting vertex and root: " << vertex.int_id() << ", " << root_vertex.int_id());
        graph_parts.vertex_to_root_.insert({vertex, root_vertex});
    }
    return graph_parts;
}
void DBGContractedGraphFactory::ProcessEdge(DBGContractedGraphFactory::contracted_dsu_t& graph_dsu,
                                            PartsBasedContractedFactory::ContractedGraphParts& parts,
                                            const std::unordered_map<VertexId, size_t>& vertex_to_id,
                                            const std::unordered_map<size_t, VertexId>& id_to_vertex,
                                            const EdgeId& edge) const {
    VertexId start = g_.EdgeStart(edge);
    VertexId end = g_.EdgeEnd(edge);
    size_t start_id = vertex_to_id.at(start);
    size_t end_id = vertex_to_id.at(end);
    size_t start_root = graph_dsu.find_set(start_id);
    size_t end_root = graph_dsu.find_set(end_id);
    VertexId start_root_vertex = id_to_vertex.at(start_root);
    VertexId end_root_vertex = id_to_vertex.at(start_root);
    if (not edge_predicate_(edge)) {
        if (start_root == end_root) {
            parts.vertex_to_capacity_.at(start_root_vertex) += g_.length(edge);
        } else {
            size_t start_capacity = parts.vertex_to_capacity_[start_root_vertex];
            size_t end_capacity = parts.vertex_to_capacity_[end_root_vertex];
            graph_dsu.unite(start_root, end_root);
            VertexId new_vertex = id_to_vertex.at(graph_dsu.find_set(start_root));
            parts.vertex_to_capacity_.at(new_vertex) = start_capacity + end_capacity + g_.length(edge);
        }
    } else {
        parts.long_edge_ends_.insert(start_root_vertex);
        parts.long_edge_ends_.insert(end_root_vertex);
        parts.long_edges_.push_back(edge);
    }
}
void SubgraphContractedGraphFactory::Construct() {
    ExtractSubgraphFromContractedGraph(other_, vertices_);
}
void SubgraphContractedGraphFactory::ExtractSubgraphFromContractedGraph(const ContractedGraph& other,
                                                                        const std::unordered_set<VertexId>& vertices) {
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
}
}