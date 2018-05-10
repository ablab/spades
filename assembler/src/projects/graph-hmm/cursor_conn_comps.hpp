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
    for (const GraphCursor &next_cursor : cursor.next()) {
      auto it = cursor2id.find(next_cursor);
      if (it != cursor2id.cend()) {
        size_t &next_id = it->second;
        dsu.unite(id, next_id);
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

// vim: set ts=2 sw=2 et :
