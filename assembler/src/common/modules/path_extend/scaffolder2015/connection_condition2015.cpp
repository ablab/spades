#include "connection_condition2015.hpp"
namespace path_extend {


map <debruijn_graph::EdgeId, double> ConnectionCondition::ConnectedWith(debruijn_graph::EdgeId e, const ScaffoldingUniqueEdgeStorage& storage) const {
    auto all_edges = this->ConnectedWith(e);
    map <debruijn_graph::EdgeId, double> res;
    for (auto edge: all_edges) {
        if (storage.IsUnique(edge.first)){
            res.insert(edge);
        }
    }
    return res;
}

PairedLibConnectionCondition::PairedLibConnectionCondition(const debruijn_graph::Graph &graph,
                             shared_ptr <PairedInfoLibrary> lib,
                             size_t lib_index,
                             size_t min_read_count) :
        graph_(graph),
        lib_(lib),
        lib_index_(lib_index),
        min_read_count_(min_read_count),
        //FIXME reconsider condition; config!
        left_dist_delta_(5 * (int) lib_->GetISMax()),
        right_dist_delta_(max(5 * (int) lib_->GetIsVar(), int(lib_->GetIS()))) {
}

size_t PairedLibConnectionCondition::GetLibIndex() const {
    return lib_index_;
}

map <debruijn_graph::EdgeId, double> PairedLibConnectionCondition::ConnectedWith(debruijn_graph::EdgeId e) const {
    set <debruijn_graph::EdgeId> all_edges;
    int e_length = (int) graph_.length(e);
    lib_->FindJumpEdges(e, all_edges,  e_length - left_dist_delta_, e_length + right_dist_delta_);

    map <debruijn_graph::EdgeId, double> result;
    for (auto edge : all_edges) {
        double w = GetWeight(e, edge);
        if (edge != e && edge != graph_.conjugate(e) &&
            math::ge(w, (double) min_read_count_)) {
            result[edge] = w;
        }
    }
    return result;
}

double PairedLibConnectionCondition::GetWeight(debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2) const {
    int e_length = (int) graph_.length(e1);
    return lib_->CountPairedInfo(e1, e2, e_length - left_dist_delta_, e_length + right_dist_delta_);
}

LongReadsLibConnectionCondition::LongReadsLibConnectionCondition(const debruijn_graph::Graph &graph,
                                    size_t lib_index,
                                    size_t min_read_count, const GraphCoverageMap& cov_map):graph_(graph), lib_index_(lib_index), min_read_count_(min_read_count), cov_map_(cov_map){}

map<debruijn_graph::EdgeId, double> LongReadsLibConnectionCondition::ConnectedWith(debruijn_graph::EdgeId ) const {
    return map <debruijn_graph::EdgeId, double>();
};

bool LongReadsLibConnectionCondition::CheckPath(BidirectionalPath *path, EdgeId e1, EdgeId e2) const {
    auto pos1 = path->FindAll(e1);
    if (pos1.size() != 1) return false;
    auto pos2 = path->FindAll(e2);
    if (pos2.size() != 1) {
        if (pos2.size() >= 2) {
            DEBUG("Something went wrong:: Edge " << graph_.int_id(e2) << "is called unique but presents in path twice! first edge " << graph_.int_id(e1) << " path ");
            path->Print();
        }
        return false;
    }
    if (pos1[0] == path->Size() - 1) return false;
    return true;
}

map<debruijn_graph::EdgeId, double> LongReadsLibConnectionCondition::ConnectedWith(debruijn_graph::EdgeId e, const ScaffoldingUniqueEdgeStorage& storage) const {
    map <debruijn_graph::EdgeId, double> res;
    auto cov_paths = cov_map_.GetCoveringPaths(e);
    DEBUG("Got cov paths " << cov_paths.size());
    for (const auto path: cov_paths) {
        auto pos1 = path->FindAll(e);
        if (pos1.size() != 1) {
            DEBUG("***not unique " << graph_.int_id(e) << " len " << graph_.length(e) << "***");
            continue;
        }
        size_t pos = pos1[0];
        pos++;
        while (pos < path->Size()){
            if (storage.IsUnique(path->At(pos))) {
                if (CheckPath(path, path->At(pos1[0]), path->At(pos))) {
                    res[path->At(pos)] += path->GetWeight();
                }
                break;
            }
            pos++;
        }
    }
    DEBUG("Before prefiltering " << res.size());
    auto iter = res.begin();
    while (iter != res.end()) {
        if (iter->second < min_read_count_){
            iter = res.erase(iter);
        } else {
            iter++;
        }
    }
    DEBUG("After prefiltering" << res.size());
    return res;
}

int LongReadsLibConnectionCondition::GetMedianGap(debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2) const {
    auto cov_paths = cov_map_.GetCoveringPaths(e1);
    std::vector<pair<int, double> > h;
    for (const auto path: cov_paths) {
        if (CheckPath(path, e1, e2)) {
            auto pos1 = path->FindAll(e1);
            auto pos2 = path->FindAll(e2);
            h.push_back(make_pair(path->LengthAt(pos1[0] + 1) - path->LengthAt(pos2[0]), path->GetWeight()));
        }
    }
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
    if (h.size() == 0) {
        WARN("filtering incorrectness");
        return 0;
    }

    return h[i].first;
}

size_t LongReadsLibConnectionCondition::GetLibIndex() const {
    return lib_index_;
}

ScaffoldGraphPairedConnectionCondition::ScaffoldGraphPairedConnectionCondition(const debruijn_graph::Graph &graph,
                                                                     const set<debruijn_graph::EdgeId>& graph_edges,
                                                                     shared_ptr <PairedInfoLibrary> lib,
                                                                     size_t lib_index,
                                                                     size_t always_add,
                                                                     size_t never_add,
                                                                     double relative_threshold):
    PairedLibConnectionCondition(graph, lib, lib_index, never_add),
    graph_edges_(graph_edges),
    always_add_(always_add),
    never_add_(never_add),
    relative_threshold_(relative_threshold) {}

map <debruijn_graph::EdgeId, double> ScaffoldGraphPairedConnectionCondition::ConnectedWith(debruijn_graph::EdgeId e) const {
    set <debruijn_graph::EdgeId> all_edges;
    int e_length = (int) graph_.length(e);
    lib_->FindJumpEdges(e, all_edges,  e_length - left_dist_delta_, e_length + right_dist_delta_);

    double max_weight = 0;
    for (auto edge : all_edges) {
        if (edge != e && edge != graph_.conjugate(e)) {
            double w = GetWeight(e, edge);
            if (graph_edges_.count(edge) > 0 && math::gr(w, max_weight))
                max_weight = w;
        }
    }
    double threshold = std::max((double) never_add_, std::min((double) always_add_, max_weight * relative_threshold_));
    map <debruijn_graph::EdgeId, double> result;
    for (auto edge : all_edges) {
        double w = GetWeight(e, edge);
        if (edge != e && edge != graph_.conjugate(e) &&
            math::ge(w, threshold)) {
            result[edge] = w;
        }
    }
    return result;
}


//TODO: We use same part of index twice, is it necessary?
int PairedLibConnectionCondition::GetMedianGap(debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2) const {
    std::vector<int> distances;
    std::vector<double> weights;
    int e_length = (int) graph_.length(e1);
    lib_->CountDistances(e1, e2, distances, weights);
    std::vector<pair<int, double> >h(distances.size());
    for (size_t i = 0; i< distances.size(); i++) {
//TODO:: we make same checks twice! That's bad
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
    if (h.size() == 0) {
        WARN("filtering incorrectness");
        return 0;
    }
    return (int) round(h[i].first - e_length);
}

AssemblyGraphConnectionCondition::AssemblyGraphConnectionCondition(const debruijn_graph::Graph &g,
                    size_t max_connection_length, const ScaffoldingUniqueEdgeStorage & unique_edges) :
        g_(g), max_connection_length_(max_connection_length), interesting_edge_set_(unique_edges.GetSet()), stored_distances_() {
}

map <debruijn_graph::EdgeId, double> AssemblyGraphConnectionCondition::ConnectedWith(debruijn_graph::EdgeId e) const {
    VERIFY_MSG(interesting_edge_set_.find(e)!= interesting_edge_set_.end(), " edge "<< e.int_id() << " not applicable for connection condition");
    if (stored_distances_.find(e) != stored_distances_.end()) {
        return stored_distances_[e];
    }
    stored_distances_.insert(make_pair(e, map<debruijn_graph::EdgeId, double>()));
    for (auto connected: g_.OutgoingEdges(g_.EdgeEnd(e))) {
        if (interesting_edge_set_.find(connected) != interesting_edge_set_.end()) {
            stored_distances_[e].insert(make_pair(connected, 1));
        }
    }
    DijkstraHelper<debruijn_graph::Graph>::BoundedDijkstra dijkstra(
            DijkstraHelper<debruijn_graph::Graph>::CreateBoundedDijkstra(g_, max_connection_length_));
    dijkstra.Run(g_.EdgeEnd(e));
    for (auto v: dijkstra.ReachedVertices()) {
        for (auto connected: g_.OutgoingEdges(v)) {
            if (interesting_edge_set_.find(connected) != interesting_edge_set_.end() && dijkstra.GetDistance(v) < max_connection_length_) {
                stored_distances_[e].insert(make_pair(connected, 1));
            }
        }
    }
    return stored_distances_[e];
}
void AssemblyGraphConnectionCondition::AddInterestingEdges(func::TypedPredicate<typename Graph::EdgeId> edge_condition) {
    for (auto e_iter = g_.ConstEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
        if (edge_condition(*e_iter))
            interesting_edge_set_.insert(*e_iter);
    }
}

size_t AssemblyGraphConnectionCondition::GetLibIndex() const {
    return (size_t) - 1;
}

int AssemblyGraphConnectionCondition::GetMedianGap (debruijn_graph::EdgeId, debruijn_graph::EdgeId) const {
    return 0;
}

}
