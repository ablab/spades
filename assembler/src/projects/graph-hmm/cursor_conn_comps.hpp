#pragma once

#include <common/adt/concurrent_dsu.hpp>
#include <unordered_map>



template <typename GraphCursor>
std::vector<std::vector<GraphCursor>> cursor_connected_components(const std::vector<GraphCursor> &cursors) {
  dsu::ConcurrentDSU dsu(cursors.size());
  std::unordered_map<GraphCursor, size_t> cursor2id;
  for (size_t i = 0; i < cursors.size(); ++i) {
    cursor2id[cursors[i]] = i;
  }

  for (const auto &kv : cursor2id) {
    const GraphCursor &cursor = kv.first;
    const size_t &id = kv.second;
    for (const auto &adj : {cursor.next(), cursor.prev()}) {
      for (const GraphCursor &adj_cursor : adj) {
        auto it = cursor2id.find(adj_cursor);
        if (it != cursor2id.cend()) {
          size_t &adj_id = it->second;
          dsu.unite(id, adj_id);
        }
      }
    }
  }

  std::vector<std::vector<size_t>> comps_ids;
  dsu.get_sets(comps_ids);
  std::vector<std::vector<GraphCursor>> result;
  result.reserve(comps_ids.size());
  for (const auto &comp_ids : comps_ids) {
    std::vector<GraphCursor> comp;
    comp.reserve(comp_ids.size());
    for (size_t i : comp_ids) {
      comp.push_back(cursors[i]);
    }
    result.emplace_back(std::move(comp));
  }

  return result;
}

template <typename T>
bool in_vector(const T &val, const std::vector<T> &vec) {
  return std::find(vec.cbegin(), vec.cend(), val) != vec.cend();
}

template <typename GraphCursor>
bool check_cursor_symmetry(const GraphCursor &cursor) {
  for (const auto &next_cursor : cursor.next()) {
    auto prevs = next_cursor.prev();
    if (!in_vector(cursor, prevs)) {
      ERROR(cursor << ", next: " << next_cursor << ", prevs: " << prevs);
      return false;
    }
  }
  for (const auto &prev_cursor : cursor.prev()) {
    auto nexts = prev_cursor.next();
    if (!in_vector(cursor, nexts)) {
      ERROR(cursor << ", prev " << prev_cursor << ", nexts" << nexts);
      return false;
    }
  }

  return true;
}

template <typename GraphCursor>
bool check_path_continuity(const std::vector<GraphCursor> &path) {
  for (size_t i = 1; i < path.size(); ++i) {
    auto nexts = path[i - 1].next();
    auto prevs = path[i].prev();
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

// vim: set ts=2 sw=2 et :
