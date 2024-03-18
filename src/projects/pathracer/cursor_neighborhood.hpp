//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2018-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <unordered_map>
#include <unordered_set>
#include "utils/logger/logger.hpp"

namespace impl {

template <typename GraphCursor>
std::vector<GraphCursor> depth_subset_dirty_heap(const std::vector<std::pair<GraphCursor, size_t>> &initial,
                                                 typename GraphCursor::Context context,
                                                 bool forward = true) {
    std::unordered_set<GraphCursor> visited;
    struct CursorWithDepth {
        GraphCursor cursor;
        size_t depth;
        bool operator<(const CursorWithDepth &other) const { return depth < other.depth; }
    };

    std::priority_queue<CursorWithDepth> q;

    std::unordered_map<GraphCursor, size_t> initial_map;
    for (const auto &cursor_with_depth : initial) {
        initial_map[cursor_with_depth.first] = std::max(initial_map[cursor_with_depth.first], cursor_with_depth.second);
    }

    for (const auto &cursor_with_depth : initial_map) {
        q.push({cursor_with_depth.first, cursor_with_depth.second});
    }

    INFO("Initial queue size: " << q.size());
    size_t step = 0;
    while (!q.empty()) {
        CursorWithDepth cursor_with_depth = q.top();
        q.pop();

        if (step % 1000000 == 0) {
            INFO("Step " << step << ", queue size: " << q.size() << " depth: " << cursor_with_depth.depth
                         << " visited size:" << visited.size());
        }
        ++step;

        if (visited.count(cursor_with_depth.cursor)) {
            continue;
        }

        visited.insert(cursor_with_depth.cursor);

        if (cursor_with_depth.depth > 0) {
            auto cursors = forward ? cursor_with_depth.cursor.next(context) : cursor_with_depth.cursor.prev(context);
            for (const auto &cursor : cursors) {
                if (!visited.count(cursor)) {
                    q.push({cursor, cursor_with_depth.depth - 1});
                }
            }
        }
    }

    return std::vector<GraphCursor>(visited.cbegin(), visited.cend());
}

}  // namespace impl

template <typename GraphCursor>
std::vector<GraphCursor> neighborhood(const GraphCursor &cursor, size_t depth, typename GraphCursor::Context context,
                                bool forward = true) {
    return impl::depth_subset_dirty_heap(std::vector<std::pair<GraphCursor, size_t>>({std::make_pair(cursor, depth)}), context,
                                         forward);
}

template <typename GraphCursor>
std::vector<GraphCursor> neighborhood(const std::vector<std::pair<GraphCursor, size_t>> &initial,
                                typename GraphCursor::Context context, bool forward = true) {
    return impl::depth_subset_dirty_heap(initial, context, forward);
}
