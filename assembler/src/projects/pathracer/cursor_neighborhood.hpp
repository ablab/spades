//***************************************************************************
//* Copyright (c) 2018-2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <boost/heap/binomial_heap.hpp>
#include <unordered_map>
#include <unordered_set>
#include "utils/logger/logger.hpp"

namespace impl {

template <class Key, class Value>
class IndexedPriorityQueue {
    struct KeyValue {
        Key key;
        Value value;
        bool operator<(const KeyValue &other) const { return value < other.value; }
    };

    using PriorityQueue = boost::heap::binomial_heap<KeyValue>;

public:
    void push_or_increase(const Key &k, const Value &v) {
        if (has(k)) {
            increase(k, v);
        } else {
            push(k, v);
        }
    }

    void push(const Key &k, const Value &v) { queue_.push({k, v}); }

    void increase(const Key &k, const Value &v) {
        auto handle = index_[k];
        KeyValue kv{k, v};
        queue_.increase(handle, kv);
    }

    bool has(const Key &k) const { return index_.count(k); }

    auto top() const { return queue_.top(); }

    void pop() { queue_.pop(); }

    bool empty() const { return queue_.empty(); }

    size_t size() const { return queue_.size(); }

private:
    PriorityQueue queue_;
    std::unordered_map<Key, typename PriorityQueue::handle_type> index_;
};

template <typename GraphCursor>
std::vector<GraphCursor> depth_subset(const std::vector<std::pair<GraphCursor, size_t>> &initial,
                                      typename GraphCursor::Context context, bool forward = true) {
    IndexedPriorityQueue<GraphCursor, size_t> q;

    std::unordered_map<GraphCursor, size_t> initial_map;
    for (const auto &cursor_with_depth : initial) {
        initial_map[cursor_with_depth.first] = std::max(initial_map[cursor_with_depth.first], cursor_with_depth.second);
    }

    for (auto &&cursor_with_depth : initial_map) {
        q.push(cursor_with_depth.first, cursor_with_depth.second);
        // q.push(std::move(cursor_with_depth.first), std::move(cursor_with_depth.second));
    }
    initial_map.clear();

    INFO("Initial queue size: " << q.size());
    std::unordered_set<GraphCursor> visited;
    while (!q.empty()) {
        auto cursor_with_depth = q.top();
        q.pop();

        if (visited.count(cursor_with_depth.key)) {
            continue;
        }
        visited.insert(cursor_with_depth.key);

        if (cursor_with_depth.value > 0) {
            auto cursors = forward ? cursor_with_depth.key.next(context) : cursor_with_depth.key.prev(context);
            for (const auto &cursor : cursors) {
                if (!visited.count(cursor)) {
                    q.push_or_increase(cursor, cursor_with_depth.value - 1);
                }
            }
        }
    }

    return std::vector<GraphCursor>(visited.cbegin(), visited.cend());
}

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
std::vector<GraphCursor> neighbourhood(const GraphCursor &cursor, size_t depth, typename GraphCursor::Context context,
                                bool forward = true) {
    return impl::depth_subset_dirty_heap(std::vector<std::pair<GraphCursor, size_t>>({std::make_pair(cursor, depth)}), context,
                                         forward);
}

template <typename GraphCursor>
std::vector<GraphCursor> neighbourhood(const std::vector<std::pair<GraphCursor, size_t>> &initial,
                                typename GraphCursor::Context context, bool forward = true) {
    return impl::depth_subset_dirty_heap(initial, context, forward);
}
