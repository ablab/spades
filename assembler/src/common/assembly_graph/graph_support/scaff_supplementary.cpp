#include "scaff_supplementary.hpp"
#include <algorithm>

using namespace std;
namespace path_extend {


void ScaffoldingUniqueEdgeAnalyzer::SetCoverageBasedCutoff() {
    vector <pair<double, size_t>> coverages;
    map <EdgeId, size_t> long_component;
    size_t total_len = 0, short_len = 0, cur_len = 0;

    for (auto iter = gp_.g.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
        if (gp_.g.length(*iter) > length_cutoff_) {
            coverages.push_back(make_pair(gp_.g.coverage(*iter), gp_.g.length(*iter)));
            total_len += gp_.g.length(*iter);
            long_component[*iter] = 0;
        } else {
            short_len += gp_.g.length(*iter);
        }
    }
    if (total_len == 0) {
        WARN("not enough edges longer than "<< length_cutoff_);
        return;
    }
    sort(coverages.begin(), coverages.end());
    size_t i = 0;
    while (cur_len < total_len / 2 && i < coverages.size()) {
        cur_len += coverages[i].second;
        i++;
    }
    median_coverage_ = coverages[i].first;
}


void ScaffoldingUniqueEdgeAnalyzer::FillUniqueEdgeStorage(ScaffoldingUniqueEdgeStorage &storage_) {
    storage_.unique_edges_.clear();
    size_t total_len = 0;
    size_t unique_len = 0;
    size_t unique_num = 0;
    storage_.SetMinLength(length_cutoff_);
    for (auto iter = gp_.g.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
        size_t tlen = gp_.g.length(*iter);
        total_len += tlen;
        if (gp_.g.length(*iter) >= length_cutoff_ && gp_.g.coverage(*iter) > median_coverage_ * (1 - relative_coverage_variation_)
                && gp_.g.coverage(*iter) < median_coverage_ * (1 + relative_coverage_variation_) ) {
            storage_.unique_edges_.insert(*iter);
            unique_len += tlen;
            unique_num ++;
        }
    }
    for (auto iter = storage_.begin(); iter != storage_.end(); ++iter) {
        DEBUG (gp_.g.int_id(*iter) << " " << gp_.g.coverage(*iter) << " " << gp_.g.length(*iter) );
    }
    INFO ("With length cutoff: " << length_cutoff_ <<", median long edge coverage: " << median_coverage_ << ", and maximal unique coverage: " <<
                                                                                                            relative_coverage_variation_);
    INFO("Unique edges quantity: " << unique_num << ", unique edges length " << unique_len <<", total edges length " << total_len);
    if (unique_len * 2 < total_len) {
        WARN("Less than half of genome in unique edges!");
    }

}

bool ScaffoldingUniqueEdgeAnalyzer::ConservativeByLength(EdgeId e) {
    return gp_.g.length(e) >= length_cutoff_;
}

map<EdgeId, size_t> ScaffoldingUniqueEdgeAnalyzer::FillNextEdgeVoting(BidirectionalPathMap<size_t>& active_paths, int direction) const {
    map<EdgeId, size_t> voting;
    for (const auto &pair: active_paths) {
        int current_pos = int(pair.second) + direction;
        auto path_iter = pair.first;
        //not found
        active_paths[path_iter] = path_iter->Size();
        while (current_pos >= 0 && current_pos < (int) path_iter->Size()) {
            if (gp_.g.length(path_iter->At(current_pos)) >= length_cutoff_) {
                voting[path_iter->At(current_pos)] += size_t(round(path_iter->GetWeight()));
                active_paths[path_iter] = size_t(current_pos);
                break;
            }
            current_pos += direction;
        }
    }
    return voting;
}

bool ScaffoldingUniqueEdgeAnalyzer::ConservativeByPaths(EdgeId e, shared_ptr<GraphCoverageMap> long_reads_cov_map, const pe_config::LongReads lr_config, int direction) const {
    BidirectionalPathSet all_set = long_reads_cov_map->GetCoveringPaths(e);
    BidirectionalPathMap<size_t> active_paths;
    size_t loop_weight = 0;
    size_t nonloop_weight = 0;
    DEBUG ("Checking " << gp_.g.int_id(e) <<" dir "<< direction );
    for (auto path_iter: all_set) {
        auto pos = path_iter->FindAll(e);
        if (pos.size() > 1)
//TODO:: path weight should be size_t?
            loop_weight += size_t(round(path_iter->GetWeight()));
        else {
            if (path_iter->Size() > 1) nonloop_weight += size_t(round(path_iter->GetWeight()));
            active_paths[path_iter] = pos[0];
        }
    }
//TODO: small plasmid, paths a-b-a, b-a-b ?
    if (loop_weight > 1) 
            return false;
        else
            DEBUG (gp_.g.int_id(e) << " loop/nonloop weight " << loop_weight << " " << nonloop_weight);
            
    EdgeId prev_unique = e;
    while (active_paths.size() > 0) {
        size_t alt = 0;
        size_t maxx = 0;
        map<EdgeId, size_t> voting = FillNextEdgeVoting(active_paths, direction);

        if (voting.size() == 0)
            break;
        EdgeId next_unique = prev_unique;
        for (const auto &pair: voting)
            if (pair.second > maxx) {
                next_unique = pair.first;
                maxx = pair.second;
            }
        for (const auto &pair: voting)
            //TODO:: 1 from config?
            if (pair.first != next_unique && pair.second > 1)
                alt += pair.second;
        if (maxx < lr_config.unique_edge_priority * double(alt)) {
            DEBUG("edge " << gp_.g.int_id(e) <<" dir "<< direction << " was not unique" );
            DEBUG("current edge " << gp_.g.int_id(next_unique));
            DEBUG("Paths " << active_paths.size());
            return false;
        } else {
            DEBUG("cur " << gp_.g.int_id(prev_unique) << " next " << gp_.g.int_id(next_unique) <<" sz " << active_paths.size());
            for (auto iter = active_paths.begin(); iter != active_paths.end();) {
                if (iter->second >= iter->first->Size() || iter->first->At(iter->second) != next_unique) {
                    iter = active_paths.erase(iter);
                } else {
                    iter++;
                }
            }
            prev_unique = next_unique;
            DEBUG(active_paths.size() << " "<< gp_.g.int_id(next_unique));
        }
    }
    DEBUG("edge " << gp_.g.int_id(e) <<" dir "<< direction << " was unique" );
    return true;
}

bool ScaffoldingUniqueEdgeAnalyzer::ConservativeByPaths(EdgeId e, shared_ptr<GraphCoverageMap> long_reads_cov_map, const pe_config::LongReads lr_config) const{
    return (ConservativeByPaths(e, long_reads_cov_map, lr_config, 1) && ConservativeByPaths(e, long_reads_cov_map, lr_config, -1));
}


void ScaffoldingUniqueEdgeAnalyzer::CheckCorrectness(ScaffoldingUniqueEdgeStorage& unique_storage_pb) {
    for (auto iter = gp_.g.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
        EdgeId e = *iter;
        bool e_unique = unique_storage_pb.IsUnique(e);
        bool e_conj_unique = unique_storage_pb.IsUnique(gp_.g.conjugate(e));
        VERIFY_MSG(!((e_unique && !e_conj_unique) || (!e_unique && e_conj_unique)), "Edge " << gp_.g.int_id(e) << " is not symmetrically unique with it conjugate");
        if (ConservativeByLength(e)) {
            if (e_unique) {
                DEBUG("edge " << gp_.g.int_id(e) << "is unique");
            } else {
                DEBUG("edge " << gp_.g.int_id(e) << "is not unique");
            }
        }
    }
}

set<VertexId> ScaffoldingUniqueEdgeAnalyzer::GetChildren(VertexId v, map <VertexId, set<VertexId>> &dijkstra_cash_) const {
    DijkstraHelper<debruijn_graph::Graph>::BoundedDijkstra dijkstra(
            DijkstraHelper<debruijn_graph::Graph>::CreateBoundedDijkstra(gp_.g, max_dijkstra_depth_, max_dijkstra_vertices_));
    dijkstra.Run(v);

    if (dijkstra_cash_.find(v) == dijkstra_cash_.end()) {
        auto tmp = dijkstra.ReachedVertices();
        tmp.push_back(v);
        dijkstra_cash_[v] = set<VertexId> (tmp.begin(), tmp.end());
    }
    return dijkstra_cash_[v];
}

bool ScaffoldingUniqueEdgeAnalyzer::FindCommonChildren(EdgeId e1, EdgeId e2, map <VertexId, set<VertexId>> &dijkstra_cash_) const {
    auto s1 = GetChildren(gp_.g.EdgeEnd(e1), dijkstra_cash_);
    auto s2 = GetChildren(gp_.g.EdgeEnd(e2), dijkstra_cash_);
    if (s1.find(gp_.g.EdgeStart(e2)) != s1.end()) {
        return true;
    }
    if (s2.find(gp_.g.EdgeStart(e1)) != s2.end()) {
        return true;
    }
    for (VertexId v: s1) {
        if (s2.find(v) != s2.end()) {
            DEBUG("bulge-like structure, edges "<< gp_.g.int_id(e1) << " " << gp_.g.int_id(e2));
            return true;
        }
    }
    return false;
}

bool ScaffoldingUniqueEdgeAnalyzer::FindCommonChildren(vector<pair<EdgeId, double>> &next_weights) const {
    map <VertexId, set<VertexId>> dijkstra_cash_;
    for (size_t i = 0; i < next_weights.size(); i ++) {
        for (size_t j = i + 1; j < next_weights.size(); j++) {
            if (next_weights[i].second * overwhelming_majority_ > next_weights[j].second
            && next_weights[j].second * overwhelming_majority_ > next_weights[i].second &&
                !FindCommonChildren(next_weights[i].first, next_weights[j].first, dijkstra_cash_)) {
                DEBUG("multiple paired info on edges " <<next_weights[i].first <<" and "<< next_weights[j].first);
                return false;
            }
        }
    }
    return true;
}

bool ScaffoldingUniqueEdgeAnalyzer::FindCommonChildren(EdgeId from, size_t lib_index) const{
    DEBUG("processing unique edge " << gp_.g.int_id(from));
    auto next_edges = gp_.clustered_indices[lib_index].Get(from);
    vector<pair<EdgeId, double>> next_weights;
    for (auto hist_pair: next_edges) {
        if (hist_pair.first == from || hist_pair.first == gp_.g.conjugate(from))
            continue;
        double total_w = 0;
        for (auto w: hist_pair.second)
            total_w += w.weight;
        if (math::gr(total_w, 1.0))
            next_weights.push_back(make_pair(hist_pair.first, total_w));
    }
    sort(next_weights.begin(), next_weights.end(), [&](pair<EdgeId, double>a, pair<EdgeId, double>b){
        return math::gr(a.second, b.second);
    });
//most popular edges. think whether it can be done faster
    if (next_weights.size() > max_different_edges_) {
        DEBUG(next_weights.size() << " continuations");
        next_weights.resize(max_different_edges_);
    }
    return FindCommonChildren(next_weights);
}


void ScaffoldingUniqueEdgeAnalyzer::ClearLongEdgesWithPairedLib(size_t lib_index, ScaffoldingUniqueEdgeStorage &storage_) const {
    set<EdgeId> to_erase;
    for (EdgeId edge: storage_ ) {
        if (!FindCommonChildren(edge, lib_index)) {
            to_erase.insert(edge);
            to_erase.insert(gp_.g.conjugate(edge));
        }
    }
    for (auto iter = storage_.begin(); iter !=  storage_.end(); ){
        if (to_erase.find(*iter) != to_erase.end()){
            iter = storage_.erase(iter);
        } else {
            iter++;
        }
    }
}


void ScaffoldingUniqueEdgeAnalyzer::FillUniqueEdgesWithLongReads(shared_ptr<GraphCoverageMap> long_reads_cov_map, ScaffoldingUniqueEdgeStorage& unique_storage_pb, const pe_config::LongReads lr_config) {
    for (auto iter = gp_.g.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
        EdgeId e = *iter;
        if (ConservativeByLength(e) && ConservativeByPaths(e, long_reads_cov_map, lr_config)) {
            unique_storage_pb.unique_edges_.insert(e);
        }
    }
    CheckCorrectness(unique_storage_pb);
}


}
