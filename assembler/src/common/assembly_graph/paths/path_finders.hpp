#pragma once

#include "assembly_graph/core/directions.hpp"

namespace omnigraph {
template<class Graph>
class UniquePathFinder {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    const Graph& graph_;
public:
    //todo use length bound if needed
    UniquePathFinder(const Graph& graph, size_t /*length_bound*/ =
    std::numeric_limits<size_t>::max())
            : graph_(graph) {}

    std::vector<EdgeId> operator()(EdgeId e,
                                   const AbstractDirection<Graph> &direction) const {
        std::vector<EdgeId> answer;
        EdgeId curr = e;
        answer.push_back(curr);
        std::set<EdgeId> was;
        while (direction.CheckUniqueOutgoingEdge(direction.EdgeEnd(curr))) {
            curr = direction.GetUniqueOutgoingEdge(direction.EdgeEnd(curr));
            if (was.count(curr) > 0)
                break;
            was.insert(curr);
            answer.push_back(curr);
        }
        return answer;
    }

    std::vector<EdgeId> UniquePathForward(EdgeId e) const {
        return this->operator()(e, ForwardDirection<Graph>(graph_));
    }

    std::vector<EdgeId> UniquePathBackward(EdgeId e) const {
        auto tmp = this->operator()(e, BackwardDirection<Graph>(graph_));
        return std::vector<EdgeId>(tmp.rbegin(), tmp.rend());
    }

};

template<class Graph>
class TrivialPathFinder {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

public:
    TrivialPathFinder(const Graph&, size_t = 0) {}

    std::vector<EdgeId> operator()(EdgeId e, const AbstractDirection<Graph> &) const {
        return {e};
    }

};

template<class Graph>
class PlausiblePathFinder {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    //todo remove graph_ field???
    const Graph& graph_;
    const size_t length_bound_;

    class DFS {
    private:
        const Graph &graph_;
        const AbstractDirection<Graph> &direction_;
        const size_t length_bound_;

        std::pair<size_t, EdgeId> find(EdgeId edge, size_t length) {
            length += graph_.length(edge);
            VertexId cross = direction_.EdgeEnd(edge);
            auto result = make_pair(length, edge);
            if (length < length_bound_
                && direction_.CheckUniqueIncomingEdge(cross)) {
                std::vector<EdgeId> outgoing = direction_.OutgoingEdges(cross);
                for (auto it = outgoing.begin(); it != outgoing.end(); ++it) {
                    auto candidate = find(*it, length);
                    if (candidate.first > result.first)
                        result = candidate;
                }
            }
            return result;
        }

        std::vector<EdgeId> RestoreAnswer(EdgeId start, EdgeId end) {
            std::vector<EdgeId> result;
            while (end != start) {
                result.push_back(end);
                end = direction_.GetUniqueIncomingEdge(direction_.EdgeStart(end));
            }
            result.push_back(start);
            return std::vector<EdgeId>(result.rbegin(), result.rend());
        }

    public:
        DFS(const Graph &graph, const AbstractDirection<Graph> &direction,
            size_t length_bound)
                : graph_(graph),
                  direction_(direction),
                  length_bound_(length_bound) {
        }

        std::vector<EdgeId> find(EdgeId edge) {
            return RestoreAnswer(edge, find(edge, 0).second);
        }
    };

public:
    PlausiblePathFinder(const Graph& graph, size_t length_bound)
            : graph_(graph),
              length_bound_(length_bound) {}

    std::vector<EdgeId> operator()(EdgeId e,
                                   const AbstractDirection<Graph> &direction) const {
        return DFS(graph_, direction, length_bound_).find(e);
    }

};
}