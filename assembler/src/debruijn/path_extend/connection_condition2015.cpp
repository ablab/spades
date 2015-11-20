

#include "connection_condition2015.hpp"
namespace path_extend {

    PairedLibConnectionCondition::PairedLibConnectionCondition(const debruijn_graph::Graph &graph,
                                 shared_ptr <PairedInfoLibrary> lib,
                                 size_t lib_index,
                                 size_t min_read_count) :
            graph_(graph),
            lib_(lib),
            lib_index_(lib_index),
            min_read_count_(min_read_count),
//TODO reconsider condition
            left_dist_delta_(5 * (int) lib_->GetIsVar()),
            right_dist_delta_(5 * (int) lib_->GetISMax()) {
    }

    size_t PairedLibConnectionCondition::GetLibIndex() const {
        return lib_index_;
    }

    set <debruijn_graph::EdgeId> PairedLibConnectionCondition::ConnectedWith(debruijn_graph::EdgeId e) const {
        set <debruijn_graph::EdgeId> all_edges;
        int e_length = (int) graph_.length(e);
        lib_->FindJumpEdges(e, all_edges, e_length - left_dist_delta_, e_length + right_dist_delta_);

        set <debruijn_graph::EdgeId> result;
        for (auto edge : all_edges) {
            if (edge != e && edge != graph_.conjugate(e) &&
                math::ge(GetWeight(e, edge), (double) min_read_count_)) {
                result.insert(edge);
            }
        }
        return result;
    }

    double PairedLibConnectionCondition::GetWeight(debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2) const {
        int e_length = (int) graph_.length(e1);
        double res = lib_->CountPairedInfo(e1, e2, e_length - left_dist_delta_, e_length + right_dist_delta_);
        VERIFY(res == lib_->CountPairedInfo(graph_.conjugate(e2), graph_.conjugate(e1),
                                            (int) graph_.length(e2) - left_dist_delta_,
                                            (int) graph_.length(e2) + right_dist_delta_));

        return res;
    }

//TODO: We use same part of index twice, is it necessary?
    int PairedLibConnectionCondition::GetMedianGap(debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2) const {
        std::vector<int> distances;
        std::vector<double> weights;
        int e_length = (int) graph_.length(e1);
        lib_->CountDistances(e1, e2, distances, weights);
        std::vector<pair<int, double> >h(distances.size());
        for (size_t i = 0; i< distances.size(); i++) {
            if (distances[i] >= e_length - left_dist_delta_ && distances[i] <= e_length + right_dist_delta_)
                h.push_back(std::make_pair(distances[i], weights[i]));
        }
//TODO: is it really necessary?
        std::sort(h.begin(), h.end());
        double sum = 0.0;
        double sum2 = 0.0;
        for (size_t j = 0; j< h.size(); ++j) {
            sum += h[j].second;
        }
        size_t i = 0;
        for (; i < h.size(); ++i) {
            sum2 += h[i].second;
            if (sum2 * 2 > sum)
                break;
        }
        if (i >= h.size()) {
            WARN("Count median error");
            i = h.size() - 1;
        }
        return (int) round(h[i].first - e_length);
    }

    AssemblyGraphConnectionCondition::AssemblyGraphConnectionCondition(const debruijn_graph::Graph &g, size_t max_connection_length) :
            g_(g),
            max_connection_length_(max_connection_length) {
    }

    set <debruijn_graph::EdgeId> AssemblyGraphConnectionCondition::ConnectedWith(debruijn_graph::EdgeId e) const {
        set <debruijn_graph::EdgeId> result;

        for (auto connected: g_.OutgoingEdges(g_.EdgeEnd(e))) {
            result.insert(connected);
        }
//TODO: optimization possible. Precompute all pairs of interesting connected vertex.
        DijkstraHelper<debruijn_graph::Graph>::BoundedDijkstra dijkstra(
                DijkstraHelper<debruijn_graph::Graph>::CreateBoundedDijkstra(g_, max_connection_length_));
        dijkstra.Run(g_.EdgeEnd(e));
        for (auto v: dijkstra.ReachedVertices()) {
            for (auto connected: g_.OutgoingEdges(v)) {
                result.insert(connected);
            }
        }

        return result;
    }

    double AssemblyGraphConnectionCondition::GetWeight(debruijn_graph::EdgeId, debruijn_graph::EdgeId) const {
        return 1.0;
    }

    size_t AssemblyGraphConnectionCondition::GetLibIndex() const {
        return (size_t) - 1;
    }

}
