#pragma once

#include <modules/algorithms/path_extend/weight_counter.hpp>
#include "barcode_mapper.hpp"

namespace tslr_resolver {

    typedef debruijn_graph::ConjugateDeBruijnGraph Graph;
    typedef Graph::EdgeId EdgeId;
    typedef Graph::VertexId VertexId;



    class BoundedBFS {
    private:
        const Graph& g_;
        VertexId start_;
        size_t length_bound_;
        size_t edge_threshold_;
        std::set <std::pair <EdgeId, int64_t> > result_;

    public:
        BoundedBFS(const Graph& g, VertexId start, size_t length_bound, size_t edge_threshold) :
                g_(g), start_(start), length_bound_(length_bound), edge_threshold_(edge_threshold)
                {}

        void run () {
            std::queue <std::pair<VertexId, size_t> > vertices;
            vertices.push(std::make_pair(start_, 0));
            while (!vertices.empty()) {
                auto current_vertex = vertices.front().first;
                size_t current_distance = vertices.front().second;
                vertices.pop();
                auto container = g_.OutgoingEdges(current_vertex);
                for (auto edge : container) {
                    if (g_.length(edge) > edge_threshold_) {
                        result_.insert(std::make_pair(edge, current_distance));
                    }
                    if (g_.length(edge) + current_distance < length_bound_) {
                        vertices.push(std::make_pair(g_.EdgeEnd(edge), current_distance + g_.length(edge)));
                    }

                }
            }
        }

        std::set <std::pair <EdgeId, int64_t> > ReturnResult () const {
            return result_;
        }
    };
}
