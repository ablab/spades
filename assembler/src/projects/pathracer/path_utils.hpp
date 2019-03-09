//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <aa_cursor.hpp>

template <typename T>
bool in_vector(const T &val, const std::vector<T> &vec) {
    return std::find(vec.cbegin(), vec.cend(), val) != vec.cend();
}

template <typename GraphCursor>
bool check_cursor_symmetry(const GraphCursor &cursor, typename GraphCursor::Context context) {
    for (const auto &next_cursor : cursor.next(context)) {
        auto prevs = next_cursor.prev(context);
        if (!in_vector(cursor, prevs)) {
            ERROR(cursor << ", next: " << next_cursor << ", prevs: " << prevs);
            return false;
        }
    }
    for (const auto &prev_cursor : cursor.prev(context)) {
        auto nexts = prev_cursor.next(context);
        if (!in_vector(cursor, nexts)) {
            ERROR(cursor << ", prev " << prev_cursor << ", nexts" << nexts);
            return false;
        }
    }

    return true;
}

template <typename GraphCursor>
bool check_path_continuity(const std::vector<GraphCursor> &path, typename GraphCursor::Context context) {
    for (size_t i = 1; i < path.size(); ++i) {
        auto nexts = path[i - 1].next(context);
        auto prevs = path[i].prev(context);
        if (!in_vector(path[i], nexts) || !in_vector(path[i - 1], prevs)) {
            return false;
        }
    }

    return true;
}

template <typename GraphCursor>
std::vector<GraphCursor> to_nucl_path(const std::vector<GraphCursor> &path) {
    return path;
}

template <typename GraphCursor>
std::vector<GraphCursor> to_nucl_path(const std::vector<AAGraphCursor<GraphCursor>> &path) {
    std::vector<GraphCursor> result;
    for (const auto aa_cursor : path) {
        for (const GraphCursor &cursor : aa_cursor.nucl_cursors()) {
            result.push_back(cursor);
        }
    }
    return result;
}

template <typename... Ts> using void_t = void;

template <typename T, typename = void>
struct has_edge_method : std::false_type {};

template <typename T>
struct has_edge_method<T, void_t<decltype(T{}.edge())>> : std::true_type {};

template<class GraphCursor>
std::enable_if_t<has_edge_method<GraphCursor>::value, std::vector<typename GraphCursor::EdgeId>> to_path(const std::vector<GraphCursor> &cpath) {
    std::vector<typename GraphCursor::EdgeId> path;

    size_t prev_position = 0;
    for (auto cursor : cpath) {
        const auto e = cursor.edge();
        size_t position = cursor.position();
        if (path.empty() || e != path.back() || prev_position >= position) {
            path.push_back(e);
        }
        prev_position = position;
    }

    return path;
}

template<class GraphCursor>
std::vector<typename GraphCursor::EdgeId> to_path(const std::vector<AAGraphCursor<GraphCursor>> &cpath) {
    return to_path(to_nucl_path(cpath));
}
