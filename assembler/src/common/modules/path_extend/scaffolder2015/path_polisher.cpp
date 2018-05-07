
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

PathContainer PathPolisher::PolishPaths(const PathContainer &paths) {
    PathContainer result;
    for (const auto& path_pair : paths) {
        BidirectionalPath path = Polish(*path_pair.first);
        BidirectionalPath *conjugate_path = new BidirectionalPath(Polish(path.Conjugate()));
        BidirectionalPath *re_path = new BidirectionalPath(conjugate_path->Conjugate());
        result.AddPair(re_path, conjugate_path);
    }
    InfoAboutGaps(result);
    return result;
}

size_t DijkstraGapCloser::MinPathLength(const PathsT& paths) const {
    size_t shortest_len = omnigraph::CumulativeLength(g_, paths.front());
    for (size_t j = 1; j < paths.size(); ++j) {
        size_t cur_len = omnigraph::CumulativeLength(g_, paths[j]);
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
        for (const auto& gap_closer : gap_closers_) {
            path = make_shared<BidirectionalPath>(gap_closer->CloseGaps(*path));
            if (path->Size() != prev_len){
                changed = true;
                prev_len = path->Size();
                break;
            }
        }
        count++;
        if (count > MAX_POLISH_ATTEMPTS) {
            INFO("Unexpected cycle while polishing path, stopping polishing " );
            path->PrintDEBUG();
            break;
        }
    }
    return *path;
}

BidirectionalPath PathGapCloser::CloseGaps(const BidirectionalPath &path) const {
    BidirectionalPath result(g_);
    if (path.Empty())
        return result;

    VERIFY(path.GapAt(0) == Gap());
    result.PushBack(path[0]);
    for (size_t i = 1; i < path.Size(); ++i) {
        if (g_.EdgeEnd(path[i - 1]) == g_.EdgeStart(path[i]) || path.GapAt(i).is_final) {
            result.PushBack(path[i], path.GapAt(i));
        } else {
            DEBUG("Gap between " << path[i - 1].int_id() << " and " << path[i].int_id() << " " << path.GapAt(i));
            auto new_gap = CloseGap(path, i, result);
            DEBUG("gap after " << new_gap);
            result.PushBack(path[i], new_gap);
        }
    }
    return result;
}

Gap DijkstraGapCloser::CloseGap(EdgeId target_edge, const Gap &orig_gap, BidirectionalPath &result) const {
    VertexId target_vertex = g_.EdgeStart(target_edge);
//TODO:: actually we do not need paths, only edges..
    omnigraph::PathStorageCallback<Graph> path_storage(g_);
    int process_res = omnigraph::ProcessPaths(g_, 0,
                            max_path_len_,
                            g_.EdgeEnd(result.Back()),
                            target_vertex,
                            path_storage);
    if (path_storage.size() == 0 || process_res != 0) {
//No paths found or path_processor error(in particular too many vertices in Dijkstra), keeping the gap
        DEBUG("PathProcessor nonzero exit code, gap left unchanged");
        return orig_gap;
    } else if (path_storage.size() > 1) {
//More than one result, using shortest result for gap length estimation
//We cannot use both common paths and bridges in one attempt;
        Gap gap = FillWithMultiplePaths(path_storage.paths(), result);
        if (gap == Gap::INVALID())
            gap = FillWithBridge(orig_gap, path_storage.paths(), result);
        return gap;
    } else {
//Closing the gap with the unique shortest result
        DEBUG("Unique path gap closing:");
        for (EdgeId e : path_storage.paths().front()) {
            DEBUG(e.int_id());
            result.PushBack(e);
        }
        return Gap(0);
    }
}

Gap PathExtenderGapCloser::CloseGap(EdgeId target_edge, const Gap &orig_gap, BidirectionalPath &result) const {
    size_t added = 0;
    VertexId target_vertex = g_.EdgeStart(target_edge);
    while (g_.EdgeEnd(result.Back()) != target_vertex) {
        bool has_grown = extender_->MakeGrowStep(result);
        if (!has_grown)
            break;
        added += g_.length(result.Back());
    }
    //FIXME think of checking for 0 in advance
    return Gap(orig_gap.gap - (int) added, {0, orig_gap.trash.current}, false);
}

Gap DijkstraGapCloser::FillWithBridge(const Gap &orig_gap,
                                      const PathsT& paths,
                                      BidirectionalPath& result) const {
    //TODO:: constant;
    auto counts = CountEdgesQuantity(paths, 300);
    DEBUG("filing gap with bridges");
    size_t path_quantity = paths.size();
    vector<EdgeId> bridges;
    for (const auto& pair: counts)
        if (pair.second == path_quantity)
            bridges.push_back(pair.first);

    if (bridges.empty()) {
        return orig_gap;
    } else {
        std::sort(bridges.begin(), bridges.end(), [&] (EdgeId e1, EdgeId e2) {
            return g_.length(e1) > g_.length(e2);});
        EdgeId bridge = bridges[0];

        VERIFY(orig_gap.gap >= 0 && orig_gap.NoTrash());
        int min_gap_before = orig_gap.gap;
        int min_gap_after = orig_gap.gap;
        for (const auto& path : paths) {
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
        DEBUG(bridge.int_id() << " " << min_gap_before);
        result.PushBack(bridge, Gap(min_gap_before, false));
        return Gap(min_gap_after, false);
    }
}

Gap DijkstraGapCloser::FillWithMultiplePaths(const PathsT& paths,
                                              BidirectionalPath& result) const {
    bool changed = false;
    auto left = LCP(paths);
    DEBUG("Filling gap with prefix")
    for (auto e : left) {
        DEBUG(e.int_id());
        result.PushBack(e);
        changed = true;
    }
    if (changed) {
        int gap = max(min_gap_,
                      int(MinPathLength(paths) - omnigraph::CumulativeLength(g_, left)));
        return Gap(gap, false);
    } else
        return Gap::INVALID();
}

std::map<EdgeId, size_t> DijkstraGapCloser::CountEdgesQuantity(const PathsT &paths, size_t length_limit) const {
    map<EdgeId, size_t> res;
    for (const auto& path: paths) {
        set<EdgeId> edge_set(path.begin(), path.end());
        for (const auto& e: edge_set) {
            if (g_.length(e) >= length_limit) {
                res[e] += 1;
            }
        }
    }
    return res;
};

size_t DijkstraGapCloser::MinPathSize(const PathsT& paths) const {
    size_t size = paths.front().size();
    for (size_t i = 1; i < paths.size(); ++i) {
        size = min(size, paths[i].size());
    }
    return size;
}

vector<EdgeId> DijkstraGapCloser::LCP(const PathsT& paths) const {
    bool all_equal = true;
    size_t index = 0;
    size_t min_size = MinPathSize(paths);

    while (index < min_size && all_equal) {
        for (size_t i = 1; i < paths.size(); ++i) {
            auto e = paths.front()[index];
            if (e != paths[i][index]) {
                all_equal = false;
                break;
            }
        }
        if (all_equal)
            ++index;
    }

    vector<EdgeId> result;
    for (size_t i = 0; i < index; ++i) {
        result.push_back(paths.front()[i]);
    }
    return result;
}

EdgeId MatePairGapCloser::FindNext(const BidirectionalPath& path,
                                   const set<EdgeId>& present_in_paths,
                                   VertexId last_v, EdgeId target_edge) const {
    auto next_edges = g_.OutgoingEdges(last_v);
    map<EdgeId, double> candidates;

    for (const auto edge: next_edges)
        if (present_in_paths.find(edge) != present_in_paths.end())
            candidates.insert(make_pair(edge, 0));

    if (candidates.size() <= 1) {
        if (candidates.size() == 0 || candidates.begin()->first == target_edge)
            return EdgeId();
        else 
            return (candidates.begin()->first);
    } else {
        int i = (int) path.Size() - 1;
        for (; i >= 0; i--) {
            if (storage_.IsUnique(path[i]))
                break;
        }
        if (i < 0) {
            return EdgeId();
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
            if (to_sort[0].second > to_sort[1].second * weight_priority && to_sort[0].first != target_edge)
                return to_sort[0].first;
            else
                return EdgeId();
        }
    }
}
//FIXME review logic
Gap MatePairGapCloser::CloseGap(EdgeId target_edge, const Gap &orig_gap, BidirectionalPath &path) const {
    VertexId target_vertex = g_.EdgeStart(target_edge);
//TODO:: condition about trash_previous - do we need it globally?
    if (orig_gap.gap <= min_gap_ || orig_gap.trash.previous > 0) {
        return orig_gap;
    } else {
        vector<EdgeId> addition;
        EdgeId last_e = path.Back();
        VertexId last_v = g_.EdgeEnd(last_e);
        DEBUG("Closing gap with mate pairs between edge " << g_.int_id(last_e)
                  << " and edge " << g_.int_id(target_edge) << " was " << orig_gap);
        omnigraph::PathStorageCallback<Graph> path_storage(g_);
        int process_res = omnigraph::ProcessPaths(g_, 0,
                                max_path_len_,
                                last_v,
                                target_vertex,
                                path_storage);
        if (process_res != 0) {
            DEBUG("PathProcessor nonzero exit code, gap left unchanged");
            return orig_gap;
        }
        set<EdgeId> present_in_paths;
        for (const auto &p: path_storage.paths())
            for (EdgeId e : p)
                present_in_paths.insert(e);

        size_t total = 0;
        while (last_e != EdgeId()) {
            last_e = FindNext(path, present_in_paths, last_v, target_edge);
            if (last_e != EdgeId()) {
                last_v = g_.EdgeEnd(last_e);
                addition.push_back(last_e);
                total += g_.length(last_e);
            }
            if (total > max_path_len_) {
                DEBUG("Closing result length too long: " << total);
                return orig_gap;
            }
        }

        int len = int(CumulativeLength(g_, addition));
        Gap gap(orig_gap.gap - len, false);
        if (gap.gap < min_gap_ && addition.size() > 0) {
            //todo constant
            if (orig_gap.gap * 2 < len) {
//inserted result significantly longer than estimated gap
                DEBUG("Filled len" << len);
            }
            if (g_.EdgeEnd(addition.back()) != target_vertex)
                gap = Gap(min_gap_, false);
            else
                gap = Gap(0, false);
        }
        for (EdgeId e : addition) {
            DEBUG(g_.int_id(e));
            path.PushBack(e);
        }
        return gap;
    }
}

}
