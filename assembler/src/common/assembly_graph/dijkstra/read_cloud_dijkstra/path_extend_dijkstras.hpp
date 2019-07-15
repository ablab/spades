#pragma once

#include <map>
#include <vector>

#include "common/modules/path_extend/extension_chooser.hpp"

namespace omnigraph {
template<class Graph, typename distance_t = size_t>
class PathBacktrackingPutChecker: public VertexPutChecker<Graph, distance_t> {
 protected:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph &g_;
    const std::map<VertexId, pair<VertexId, EdgeId>> *prev_vert_map_;

 public:
    explicit PathBacktrackingPutChecker(const Graph &g) :
        VertexPutChecker<Graph, distance_t>(), g_(g), prev_vert_map_(nullptr) {}

    void AddPrevVertexMap(const std::map<VertexId, pair<VertexId, EdgeId>> &prev_vert_map) {
        prev_vert_map_ = &prev_vert_map;
    }
};

template<class Graph, typename distance_t = size_t>
class ExtensionChooserBasedPutChecker: public PathBacktrackingPutChecker<Graph, distance_t> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    using PathBacktrackingPutChecker<Graph, distance_t>::g_;
    using PathBacktrackingPutChecker<Graph, distance_t>::prev_vert_map_;

    typedef path_extend::BidirectionalPath BidirectionalPath;

    const shared_ptr<path_extend::ExtensionChooser> extension_chooser_;
    const size_t path_prefix_length_;
    const VertexId init_vertex_;
    const EdgeId init_edge_;
 public:
    ExtensionChooserBasedPutChecker(const Graph &g,
                                    const shared_ptr<path_extend::ExtensionChooser> extension_chooser_,
                                    const size_t path_prefix_length_,
                                    const VertexId &init_vertex_,
                                    const EdgeId &init_edge_)
        : PathBacktrackingPutChecker<Graph, distance_t>(g),
          extension_chooser_(extension_chooser_),
          path_prefix_length_(path_prefix_length_),
          init_vertex_(init_vertex_),
          init_edge_(init_edge_) {}

    bool Check(VertexId vertex, EdgeId edge, distance_t distance) const override {
        path_extend::ExtensionChooser::EdgeContainer outgoing_container;
        //fixme ineffective, dangerous, does not fix gaps.
        for (const EdgeId& next: g_.OutgoingEdges(g_.EdgeStart(edge))) {
            outgoing_container.emplace_back(next, 0);
        }

        DEBUG("Checking vertex " << vertex.int_id());
        DEBUG("Checking edge " << edge.int_id());
        DEBUG("Distance: " << distance);

        VertexId prev_vertex = g_.EdgeStart(edge);

        auto path_suffix = GetPathSuffix(path_prefix_length_, prev_vertex);
        DEBUG(path_suffix.Size() << " edges in preceding path");
        DEBUG("Preceding path length: " << path_suffix.Length());
        auto filtered_edges = extension_chooser_->Filter(path_suffix, outgoing_container);
        DEBUG(filtered_edges.size() << " filtered edges");

        if (filtered_edges.empty()) {
            return true;
        }

        for (const auto& passed: filtered_edges) {
            DEBUG("Comparing with candidate: " << passed.e_.int_id());
            if (edge == passed.e_) {
                DEBUG("Edge passed");
                return true;
            }
        }
        DEBUG("Edge did not pass");
        return false;
    }

    path_extend::BidirectionalPath GetPathSuffix(size_t suffix_length, VertexId vertex) const {
        DEBUG("Getting previous path from vertex: " << vertex.int_id());
        DEBUG("Initial vertex: " << init_vertex_.int_id());
        if (vertex == init_vertex_) {
            std::vector<EdgeId> init_path({init_edge_});
            BidirectionalPath result(g_, init_path);
            return result;
        }

        std::vector<EdgeId> path;
        if (prev_vert_map_->find(vertex) == prev_vert_map_->end()) {
            BidirectionalPath result(g_, path);
            return result;
        }

        VertexId curr_vertex = vertex;
        VertexId prev_vertex = utils::get(*prev_vert_map_, curr_vertex).first;
        EdgeId prev_edge = utils::get(*prev_vert_map_, curr_vertex).second;
        size_t current_suff_len = 0;

        while (prev_vertex != VertexId(0) and current_suff_len < suffix_length) {
            DEBUG("Previous vertex: " << prev_vertex.int_id());
            VERIFY(prev_vertex == g_.EdgeStart(prev_edge));
            path.push_back(prev_edge);
            curr_vertex = prev_vertex;
            DEBUG("Searching in map");
            DEBUG("Map size: " << prev_vert_map_->size());
            if (prev_vert_map_->find(curr_vertex) == prev_vert_map_->end()) {
                break;
            }
            const auto& prev_v_e = utils::get(*prev_vert_map_, curr_vertex);
            prev_vertex = prev_v_e.first;
            DEBUG("Prev vertex id: " << prev_vertex.int_id());
            prev_edge = prev_v_e.second;
            DEBUG("Prev edge id: " << prev_edge.int_id());
            if (prev_vertex.int_id() == 0) {
                if (curr_vertex == init_vertex_) {
                    path.push_back(init_edge_);
                }
                break;
            }
            current_suff_len += g_.length(prev_edge);
            DEBUG("Current suffix length: " << current_suff_len);
        }
        std::reverse(path.begin(), path.end());
        DEBUG("Suffix path: ");
        for (const auto& edge: path) {
            DEBUG(edge.int_id());
        }
        BidirectionalPath result(g_, path);
        return result;
    }

    DECL_LOGGER("ExtensionChooserBasedPutChecker");
};
}