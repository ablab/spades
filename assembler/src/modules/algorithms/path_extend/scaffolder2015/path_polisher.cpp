//
// Created by andrey on 12.05.16.
//

#include "path_polisher.hpp"

namespace path_extend {


void PathGapCloser::PolishPaths(const PathContainer &paths, PathContainer &result) {
    result.clear();

    for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
        BidirectionalPath *path = new BidirectionalPath(Polish(*iter.get()));
        BidirectionalPath *conjugatePath = new BidirectionalPath(path->Conjugate());
        result.AddPair(path, conjugatePath);
    }
}

size_t DijkstraGapCloser::MinPathLength(const omnigraph::PathStorageCallback<Graph>& path_storage) const {
    size_t shortest_len = omnigraph::CumulativeLength(g_, path_storage.paths().front());
    for (size_t j = 1; j < path_storage.paths().size(); ++j) {
        size_t cur_len = omnigraph::CumulativeLength(g_, path_storage.paths()[j]);
        shortest_len = min(shortest_len, cur_len);
    }
    return shortest_len;
}


BidirectionalPath DijkstraGapCloser::Polish(const BidirectionalPath &path) {
    BidirectionalPath result(g_);
    if (path.Empty())
        return result;

    result.PushBack(path.Front(), path.GapAt(0));
    for (size_t i = 1; i < path.Size(); ++i) {
        if (g_.EdgeEnd(path[i - 1]) == g_.EdgeStart(path[i])) {
            result.PushBack(path[i], path.GapAt(i));
        }
        else {
            //Connect edges using Dijkstra
            omnigraph::PathStorageCallback<Graph> path_storage(g_);
            omnigraph::ProcessPaths(g_, 0,
                                    max_path_len_,
                                    g_.EdgeEnd(path[i - 1]),
                                    g_.EdgeStart(path[i]),
                                    path_storage);

            if (path_storage.size() == 0) {
                //No paths found, keeping the gap
                result.PushBack(path[i], path.GapAt(i));
            }
            else if (path_storage.size() > 1) {
                //More than one path, using shortest path for gap length estimation
                FillWithMultiplePaths(path, i, path_storage, result);
            }
            else {
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

void DijkstraGapCloser::FillWithMultiplePaths(const BidirectionalPath& path, size_t index,
                                              const omnigraph::PathStorageCallback<Graph>& path_storage,
                                              BidirectionalPath& result) {

    result.PushBack(path[index], max((int) MinPathLength(path_storage), 100));
}


void CommonPrefixDijkstraGapCloser::FillWithMultiplePaths(const BidirectionalPath& path, size_t index,
                                              const omnigraph::PathStorageCallback<Graph>& path_storage,
                                              BidirectionalPath& result) {

    auto left = LCP(path_storage);
    auto right = LCS(path_storage);
    if (left.size() + right.size() > MinPathSize(path_storage)) {
        if (!right.empty()) {
            auto last_edge = right.front();
            while (!left.empty()) {
                if (left.back() == last_edge) {
                    left.pop_back();
                    break;
                }
                left.pop_back();
            }
        }
    }

    for (auto e : left) {
        result.PushBack(e);
    }

    int middle_gap = (int) max(size_t(100), MinPathLength(path_storage) -
            omnigraph::CumulativeLength(g_, left) - omnigraph::CumulativeLength(g_, right));

    if (right.empty()) {
        result.PushBack(path[index], middle_gap);
    }
    else {
        result.PushBack(right.front(), middle_gap);
        for (size_t i = 1; i < right.size(); ++i) {
            result.PushBack(right[i]);
        }
        result.PushBack(path[index]);
    }
}


size_t CommonPrefixDijkstraGapCloser::MinPathSize(const omnigraph::PathStorageCallback<Graph>& path_storage) const {
    size_t size = path_storage.paths().front().size();
    for (size_t i = 1; i < path_storage.size(); ++i) {
        size = min(size, path_storage.paths()[i].size());
    }
    return size;
}

vector<EdgeId> CommonPrefixDijkstraGapCloser::LCP(const omnigraph::PathStorageCallback<Graph>& path_storage) const {
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

vector<EdgeId> CommonPrefixDijkstraGapCloser::LCS(const omnigraph::PathStorageCallback<Graph>& path_storage) const {
    bool all_equal = true;
    size_t index = 0;
    size_t min_size = MinPathSize(path_storage);

    while (index < min_size && all_equal) {
        for (size_t i = 1; i < path_storage.size(); ++i) {
            auto e = path_storage.paths().front()[path_storage.paths().front().size() - index - 1];
            if (e != path_storage.paths()[i][path_storage.paths()[i].size() - index - 1]) {
                all_equal = false;
                break;
            }
        }
        if (all_equal)
            ++index;
    }

    vector<EdgeId> result;
    auto path = path_storage.paths().front();
    for (size_t i = path.size() - index; i < path.size(); ++i) {
        result.push_back(path[i]);
    }
    return result;
}
}