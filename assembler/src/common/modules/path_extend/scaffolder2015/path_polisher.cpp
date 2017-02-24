
#include "path_polisher.hpp"

namespace path_extend {

void PathPolisher::InfoAboutGaps(const PathContainer & result){
    for (const auto& p_iter: result) {
        for (size_t i = 1; i < p_iter.first->Size(); ++i) {
            if (p_iter.first->GapAt(i) > 0) {
                DEBUG("Gap "<< p_iter.first->GapAt(i) << " left between " << gp_.g.int_id(p_iter.first->At(i-1)) << " and " << gp_.g.int_id(p_iter.first->At(i)));
            }
        }
    }
}

PathPolisher::PathPolisher(const conj_graph_pack& gp, const config::dataset& dataset_info, const ScaffoldingUniqueEdgeStorage& storage, size_t max_resolvable_len ): gp_(gp) {
    gap_closers.push_back(make_shared<DijkstraGapCloser>(gp.g, max_resolvable_len));
    for (size_t i = 0; i <  dataset_info.reads.lib_count(); i++) {
        auto lib = dataset_info.reads[i];
        if (lib.type() == io::LibraryType::HQMatePairs || lib.type() == io::LibraryType::MatePairs) {
            shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp.g, lib, gp.paired_indices[i]);
            gap_closers.push_back(make_shared<MatePairGapCloser> (gp.g, max_resolvable_len, paired_lib, storage));
        }
    }
}

void PathPolisher::PolishPaths(const PathContainer &paths, PathContainer &result) {
    result.clear();

    for (auto iter = paths.begin(); iter != paths.end(); ++iter) {

        BidirectionalPath *path = new BidirectionalPath(Polish(*iter.get()));
        BidirectionalPath *conjugatePath = new BidirectionalPath(Polish(path->Conjugate()));
        BidirectionalPath *re_path = new BidirectionalPath(conjugatePath->Conjugate());
        result.AddPair(re_path, conjugatePath);
    }
    InfoAboutGaps(result);
}

size_t DijkstraGapCloser::MinPathLength(const omnigraph::PathStorageCallback<Graph>& path_storage) const {
    size_t shortest_len = omnigraph::CumulativeLength(g_, path_storage.paths().front());
    for (size_t j = 1; j < path_storage.paths().size(); ++j) {
        size_t cur_len = omnigraph::CumulativeLength(g_, path_storage.paths()[j]);
        shortest_len = min(shortest_len, cur_len);
    }
    return shortest_len;
}

BidirectionalPath PathPolisher::Polish(const BidirectionalPath &path) {
    if (path.Empty())
        return path;
    shared_ptr<BidirectionalPath> current;
    shared_ptr<BidirectionalPath> prev_step = std::make_shared<BidirectionalPath>(path);
    bool changed = true;
    size_t count = 0;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < gap_closers.size(); i++) {
            current = std::make_shared<BidirectionalPath>(gap_closers[i]->Polish(*prev_step));
            if (current->Size() != prev_step->Size()){
                changed = true;
                std::swap(current, prev_step);
                break;
            }
        }
        count++;
        if (count > 5) {
            INFO("Unexpected cycle while polishing path, stopping polishing " );
            path.Print();
            break;
        }
    }
    return *prev_step;
}

BidirectionalPath DijkstraGapCloser::Polish(const BidirectionalPath &path) {
    BidirectionalPath result(g_);
    if (path.Empty())
        return result;
    result.PushBack(path[0], path.GapInfoAt(0));
    for (size_t i = 1; i < path.Size(); ++i) {
        if (g_.EdgeEnd(path[i - 1]) == g_.EdgeStart(path[i])) {
            result.PushBack(path[i], path.GapInfoAt(i));
        } else {
            //Connect edges using Dijkstra
            omnigraph::PathStorageCallback<Graph> path_storage(g_);
            omnigraph::ProcessPaths(g_, 0,
                                    max_path_len_,
                                    g_.EdgeEnd(path[i - 1]),
                                    g_.EdgeStart(path[i]),
                                    path_storage);

            if (path_storage.size() == 0) {
                //No paths found, keeping the gap
                result.PushBack(path[i], path.GapInfoAt(i));
            } else if (path_storage.size() > 1) {
                //More than one path, using shortest path for gap length estimation
                //We cannot use both common paths and bridges in one attempt;
                if (!FillWithMultiplePaths(path, i, path_storage, result))
                    FillWithBridge(path, i, path_storage, result);
            } else {
                //Closing the gap with the unique shortest path
                for (size_t j = 0; j < path_storage.paths().front().size(); ++j) {
                    result.PushBack(path_storage.paths().front()[j]);
                }
                result.PushBack(path[i]);
            }
        }
    }
    return result;
}


bool DijkstraGapCloser::FillWithBridge(const BidirectionalPath& path, size_t index,
                                                          const omnigraph::PathStorageCallback<Graph>& path_storage,
                                                          BidirectionalPath& result) const {
//TODO:: constant;
    auto counts = CountEdgesQuantity(path_storage, 300);
    size_t path_quantity = path_storage.paths().size();
    vector<EdgeId> bridges;
    for (const auto& pair: counts)
        if (pair.second == path_quantity)
            bridges.push_back(pair.first);
    if (bridges.size() > 0) {
        std::sort(bridges.begin(), bridges.end(), [&] (EdgeId e1, EdgeId e2) {
            return g_.length(e1) > g_.length(e2); });
        EdgeId bridge = bridges[0];
        int min_gap_before = path.GapAt(index);
        int min_gap_after = path.GapAt(index);
        for (const auto& path:path_storage.paths()) {
            int current_before = 0;
            for(size_t i = 0; i< path.size(); i++) {
                if (path[i] != bridge)
                    current_before += (int)g_.length(path[i]);
                else
                    break;
            }
            int current_after = (int)CumulativeLength(g_, path) - current_before - int(g_.length(bridge));
            min_gap_after = std::min(current_after, min_gap_after);
            min_gap_before = std::min(current_before, min_gap_before);
        }
        min_gap_after = std::max(min_gap_after, min_gap_);
        min_gap_before = std::max(min_gap_before, min_gap_);
        result.PushBack(bridge, min_gap_before);
        result.PushBack(path[index], min_gap_after);
        return true;
    } else {
        result.PushBack(path[index], path.GapAt(index));
        return false;
    }
}

bool DijkstraGapCloser::FillWithMultiplePaths(const BidirectionalPath& path, size_t index,
                                              const omnigraph::PathStorageCallback<Graph>& path_storage,
                                              BidirectionalPath& result) const {
    bool changed = false;
    auto left = LCP(path_storage);
    for (auto e : left) {
        result.PushBack(e);
        changed = true;
    }
    int middle_gap = (int) max(size_t(min_gap_), MinPathLength(path_storage) -
            omnigraph::CumulativeLength(g_, left));
    if (changed)
        result.PushBack(path[index], middle_gap);
    return changed;
}

std::map<EdgeId, size_t> DijkstraGapCloser::CountEdgesQuantity(const omnigraph::PathStorageCallback<Graph>& path_storage, size_t length_limit ) const{
    map<EdgeId, size_t> res;
    for (const auto& path: path_storage.paths()) {
        set<EdgeId> edge_set(path.begin(), path.end());
        for (const auto& e: edge_set) {
            if (g_.length(e) >= length_limit) {
                res[e] += 1;
            }
        }
    }
    return res;
};

size_t DijkstraGapCloser::MinPathSize(const omnigraph::PathStorageCallback<Graph>& path_storage) const {
    size_t size = path_storage.paths().front().size();
    for (size_t i = 1; i < path_storage.size(); ++i) {
        size = min(size, path_storage.paths()[i].size());
    }
    return size;
}

vector<EdgeId> DijkstraGapCloser::LCP(const omnigraph::PathStorageCallback<Graph>& path_storage) const {
    bool all_equal = true;
    size_t index = 0;
    size_t min_size = MinPathSize(path_storage);

    while (index < min_size && all_equal) {
        for (size_t i = 1; i < path_storage.size(); ++i) {
            auto e = path_storage.paths().front()[index];
            if (e != path_storage.paths()[i][index]) {
                all_equal = false;
                break;
            }
        }
        if (all_equal)
            ++index;
    }

    vector<EdgeId> result;
    for (size_t i = 0; i < index; ++i) {
        result.push_back(path_storage.paths().front()[i]);
    }
    return result;
}


EdgeId MatePairGapCloser::FindNext(const BidirectionalPath& path, size_t index,
                    const set<EdgeId>& present_in_paths, VertexId v) const {
    auto next_edges = g_.OutgoingEdges(v);
    map<EdgeId, double> candidates;
    for (const auto edge: next_edges)
        if (present_in_paths.find(edge) != present_in_paths.end())
            candidates.insert(make_pair(edge, 0));
    if (candidates.size() <= 1 ) {
        if (candidates.size() == 0 || candidates.begin()->first == path[index])
            return EdgeId(0);
        else 
            return (candidates.begin()->first);
    } else {
        int i = (int) index - 1;
        for (; i >= 0; i--) {
            if (storage_.IsUnique(path[i]))
                break;
        }
        if (i < 0) {
            return EdgeId(0);
        } else {
            EdgeId last_unique = path[i];
            for (auto &pair: candidates){
                vector<int> d;
                vector<double> w;
//TODO:: any filtration?
                lib_->CountDistances(last_unique, pair.first, d, w);
                double sum = 0;
                for (auto weight: w)
                    sum += weight;
                pair.second = sum / double(g_.length(pair.first));
            }
            vector<std::pair<EdgeId, double>> to_sort(candidates.begin(),candidates.end());
            sort(to_sort.begin(), to_sort.end(), [&] (std::pair<EdgeId, double> a, std::pair<EdgeId, double> b ) {
                return a.second > b.second;
            });
            if (to_sort[0].second > to_sort[1].second * weight_priority && to_sort[0].first != path[index])
                return to_sort[0].first;
            else
                return EdgeId(0);
        }
    }
}

//TODO: make shorter functions
BidirectionalPath MatePairGapCloser::Polish(const BidirectionalPath& path) {
    BidirectionalPath result(g_);
    DEBUG("Path " << path.GetId() << " len "<< path.Length() << " size " << path.Size());
    result.PushBack(path[0], path.GapInfoAt(0));
    for (size_t i = 1; i < path.Size(); ++i) {
        if (g_.EdgeEnd(path[i - 1]) == g_.EdgeStart(path[i]) || path.GapAt(i) <= min_gap_) {
            result.PushBack(path[i], path.GapInfoAt(i));
        } else {
            DEBUG("position "<< i <<" gap between edges " << g_.int_id(path[i-1]) << " and " << g_.int_id(path[i]) << " was " << path.GapAt(i));

            vector<EdgeId> addition;
            VertexId v = g_.EdgeEnd(path[i - 1]);
            EdgeId last = path[i - 1];
            omnigraph::PathStorageCallback<Graph> path_storage(g_);
            omnigraph::ProcessPaths(g_, 0,
                                    max_path_len_,
                                    g_.EdgeEnd(path[i - 1]),
                                    g_.EdgeStart(path[i]),
                                    path_storage);
            set<EdgeId> present_in_paths;
            for(const auto &p: path_storage.paths())
                for(size_t j = 0; j < p.size(); j ++)
                    present_in_paths.insert(p[j]);
            size_t total = 0;
            while (last != EdgeId(0)){
                last = FindNext(path, i, present_in_paths, v);
                if (last != EdgeId(0)){
                    v = g_.EdgeEnd(last);
                    addition.push_back(last);
                    total += g_.length(last);
                }
                if (total > max_path_len_){
                    DEBUG("gap between edges " << g_.int_id(path[i-1]) << " and " << g_.int_id(path[i]) << " was: " << path.GapAt(i) << ", closing path length too long: " << total);
                    break;
                }
            }
            if (total > max_path_len_) {
                result.PushBack(path[i], path.GapInfoAt(i));
                continue;                
            }
            int len = int(CumulativeLength(g_, addition));
            int new_gap = path.GapAt(i) - len;
            if (new_gap < min_gap_ && addition.size() > 0) {
                if (path.GapAt(i) * 3 < len * 2 ) {
//inserted path significantly longer than estimated gap
                    DEBUG("Gap size estimation problem: gap between edges " << g_.int_id(path[i - 1]) << " and " << g_.int_id(path[i]) << " was " <<
                         path.GapAt(i) << "filled len" << len);
                }
                if (g_.EdgeEnd(addition.back()) != g_.EdgeStart(path[i]))
                    new_gap = min_gap_;
                else
                    new_gap = 0;
            }
            DEBUG("filling");
            for (size_t j = 0; j < addition.size(); j++) {
                DEBUG(g_.int_id(addition[j]));
                result.PushBack(addition[j], 0);
            }
            result.PushBack(path[i], new_gap);
            DEBUG("filled");
        }
    }
    DEBUG("result " << result.GetId() << " len "<< result.Length() << " size " << result.Size());
    return result;
}

}
