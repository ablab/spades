
#include "path_polisher.hpp"

namespace path_extend {

void PathPolisher::InfoAboutGaps(const PathContainer & result){
    for (const auto& p_iter: result) {
        for (size_t i = 1; i < p_iter.first->Size(); ++i) {
            if (p_iter.first->GapAt(i).gap > 0) {
                DEBUG("Gap "<< p_iter.first->GapAt(i).gap
                            << " left between " << gp_.g.int_id(p_iter.first->At(i-1))
                            << " and " << gp_.g.int_id(p_iter.first->At(i)));
            }
        }
    }
}

PathPolisher::PathPolisher(const conj_graph_pack& gp, const config::dataset& dataset_info,
                           const ScaffoldingUniqueEdgeStorage& storage, size_t max_resolvable_len, vector<shared_ptr<PathExtender>> extenders): gp_(gp) {
    gap_closers.push_back(make_shared<DijkstraGapCloser>(gp.g, max_resolvable_len));
    for (size_t i = 0; i <  dataset_info.reads.lib_count(); i++) {
        auto lib = dataset_info.reads[i];
        if (lib.type() == io::LibraryType::HQMatePairs || lib.type() == io::LibraryType::MatePairs) {
            shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp.g, lib, gp.paired_indices[i]);
            gap_closers.push_back(make_shared<MatePairGapCloser> (gp.g, max_resolvable_len, paired_lib, storage));
        }
    }
    for (auto extender: extenders) {
        PathExtenderGapCloser(gp.g, max_resolvable_len, extender);
    }
}

void PathPolisher::PolishPaths(const PathContainer &paths, PathContainer &result) {
    result.clear();

    for (const auto& path_pair : paths) {
        BidirectionalPath path = Polish(*path_pair.first);
        BidirectionalPath *conjugate_path = new BidirectionalPath(Polish(path.Conjugate()));
        BidirectionalPath *re_path = new BidirectionalPath(conjugate_path->Conjugate());
        result.AddPair(re_path, conjugate_path);
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

BidirectionalPath PathPolisher::Polish(const BidirectionalPath &init_path) {
    if (init_path.Empty())
        return init_path;

    auto path = make_shared<BidirectionalPath>(init_path);
    size_t prev_len = path->Size();

    bool changed = true;
    size_t count = 0;
    while (changed) {
        changed = false;
        for (auto& gap_closer : gap_closers) {
            path = make_shared<BidirectionalPath>(gap_closer->Polish(*path));
            if (path->Size() != prev_len){
                changed = true;
                prev_len = path->Size();
                break;
            }
        }
        count++;
        if (count > max_polish_attempts) {
            INFO("Unexpected cycle while polishing path, stopping polishing " );
            path->Print();
            break;
        }
    }
    return *path;
}

BidirectionalPath PathGapCloser::Polish(const BidirectionalPath &path) {
    BidirectionalPath result(g_);
    if (path.Empty())
        return result;

    VERIFY(path.GapAt(0) == Gap());
    result.PushBack(path[0]);
    for (size_t i = 1; i < path.Size(); ++i) {
        if (g_.EdgeEnd(path[i - 1]) == g_.EdgeStart(path[i])) {
            result.PushBack(path[i], path.GapAt(i));
        } else {
            auto new_gap = InnerCloseGap(path, i, result);
            result.PushBack(path[i], new_gap);
        }
    }
    return result;
}

Gap DijkstraGapCloser::InnerCloseGap(const BidirectionalPath &original_path, size_t position, BidirectionalPath &path) {
    omnigraph::PathStorageCallback<Graph> path_storage(g_);
    omnigraph::ProcessPaths(g_, 0,
                            max_path_len_,
                            g_.EdgeEnd(original_path[position - 1]),
                            g_.EdgeStart(original_path[position]),
                            path_storage);
    if (path_storage.size() == 0) {
//No paths found, keeping the gap
        return path.GapAt(position);
    } else if (path_storage.size() > 1) {
//More than one path, using shortest path for gap length estimation
//We cannot use both common paths and bridges in one attempt;
        Gap gap = FillWithMultiplePaths(path_storage, path);
        if (gap == Gap::INVALID())
            gap = FillWithBridge(original_path, position, path_storage, path);
        return gap;
    } else {
//Closing the gap with the unique shortest path
        for (EdgeId e : path_storage.paths().front()) {
            path.PushBack(e);
        }
        return Gap(0);
    }
}

Gap PathExtenderGapCloser::InnerCloseGap(const BidirectionalPath &original_path, size_t position, BidirectionalPath &path) {
    size_t added = 0;
    VertexId stop_vertex = g_.EdgeStart(original_path.At(position));
    while (g_.EdgeEnd(path.Back()) != stop_vertex) {
        bool has_grown = extender_->MakeGrowStep(path);
        if (!has_grown)
            break;
        added += g_.length(path.Back());
    }
    return Gap(path.GapAt(position).gap - (int) added, 0, path.GapAt(position).trash_current);
}



Gap DijkstraGapCloser::FillWithBridge(const BidirectionalPath& path, size_t index,
                                                          const omnigraph::PathStorageCallback<Graph>& path_storage,
                                                          BidirectionalPath& result) const {
    //TODO:: constant;
    auto counts = CountEdgesQuantity(path_storage, 300);
    size_t path_quantity = path_storage.paths().size();
    vector<EdgeId> bridges;
    for (const auto& pair: counts)
        if (pair.second == path_quantity)
            bridges.push_back(pair.first);

    if (bridges.empty()) {
        return path.GapAt(index);
    } else {
        std::sort(bridges.begin(), bridges.end(), [&] (EdgeId e1, EdgeId e2) {
            return g_.length(e1) > g_.length(e2); });
        EdgeId bridge = bridges[0];

        Gap orig_gap = path.GapAt(index);
        VERIFY(path.GapAt(index).gap >= 0 && orig_gap.NoTrash());
        int min_gap_before = path.GapAt(index).gap;
        int min_gap_after = path.GapAt(index).gap;
        for (const auto& path : path_storage.paths()) {
            size_t current_before = 0;
            for (EdgeId e : path) {
                if (e == bridge)
                    break;
                current_before += g_.length(e);
            }
            size_t current_after = CumulativeLength(g_, path) - current_before - g_.length(bridge);
            min_gap_after = std::min(int(current_after), min_gap_after);
            min_gap_before = std::min(int(current_before), min_gap_before);
        }

        min_gap_after = std::max(min_gap_after, min_gap_);
        min_gap_before = std::max(min_gap_before, min_gap_);
        result.PushBack(bridge, Gap(min_gap_before));
        return Gap(min_gap_after);
    }
}

Gap DijkstraGapCloser::FillWithMultiplePaths(const omnigraph::PathStorageCallback<Graph>& path_storage,
                                              BidirectionalPath& result) const {
    bool changed = false;
    auto left = LCP(path_storage);
    for (auto e : left) {
        result.PushBack(e);
        changed = true;
    }
    if (changed) {
        int gap = max(min_gap_, int(MinPathLength(path_storage) - omnigraph::CumulativeLength(g_, left)));
        return Gap(gap);
    } else
        return Gap::INVALID();
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

Gap MatePairGapCloser::InnerCloseGap(const BidirectionalPath &original_path, size_t position,
                                                   BidirectionalPath &path)  {
//TODO:: condition about trash_previous - do we need it globally?
    if (original_path.GapAt(position).gap <= min_gap_ || original_path.GapAt(position).trash_previous > 0) {
        return path.GapAt(position);
    } else {
        DEBUG("position "<< position <<" gap between edges " << g_.int_id(original_path.At(position - 1))
              << " and " << g_.int_id(original_path.At(position)) << " was " << original_path.GapAt(position).gap);
        vector<EdgeId> addition;
        EdgeId last = original_path.At(position - 1);
        VertexId v = g_.EdgeEnd(last);
        omnigraph::PathStorageCallback<Graph> path_storage(g_);
        omnigraph::ProcessPaths(g_, 0,
                                max_path_len_,
                                g_.EdgeEnd(original_path.At(position - 1)),
                                g_.EdgeStart(original_path.At(position)),
                                path_storage);
        set<EdgeId> present_in_paths;
        for (const auto &p: path_storage.paths())
            for (EdgeId e : p)
                present_in_paths.insert(e);

        size_t total = 0;
        while (last != EdgeId(0)) {
            last = FindNext(original_path, position, present_in_paths, v);
            if (last != EdgeId(0)){
                v = g_.EdgeEnd(last);
                addition.push_back(last);
                total += g_.length(last);
            }
            if (total > max_path_len_){
                DEBUG("gap between edges " << g_.int_id(original_path.At(position -1)) << " and " << g_.int_id(original_path.At(position))
                      << " was: " << original_path.GapAt(position).gap << ", closing path length too long: " << total);
                break;
            }
        }
        if (total > max_path_len_) {
            return path.GapAt(position);
        }
        int len = int(CumulativeLength(g_, addition));
        Gap gap(original_path.GapAt(position).gap - len);
        if (gap.gap < min_gap_ && addition.size() > 0) {
            if (original_path.GapAt(position).gap * 2 < len  ) {
//inserted path significantly longer than estimated gap
                DEBUG("Gap size estimation problem: gap between edges " << g_.int_id(original_path.At(position - 1))
                      << " and " << g_.int_id(original_path.At(position)) << " was "
                      << original_path.GapAt(position).gap << "filled len" << len);
            }
            if (g_.EdgeEnd(addition.back()) != g_.EdgeStart(original_path.At(position)))
                gap = Gap(min_gap_);
            else
                gap = Gap();
        }
        DEBUG("filling");
        for (EdgeId e : addition) {
            DEBUG(g_.int_id(e));
            path.PushBack(e);
        }
        return gap;
    }
}

}
