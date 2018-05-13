#pragma once

#include <vector>
#include <utility>
#include <unordered_set>

template <class GraphCursor>
class ReversalGraphCursor : public GraphCursor {
 public:
  using GraphCursor::GraphCursor;
  using GraphCursor::operator=;
  // "an inherited constructor is not a candidate for initialization from an expression of the same or derived type"
  ReversalGraphCursor(const GraphCursor &other) : GraphCursor(other) {}
  ReversalGraphCursor(GraphCursor &&other) : GraphCursor(std::move(other)) {}

  ReversalGraphCursor() = default;
  ReversalGraphCursor(const ReversalGraphCursor&) = default;
  ReversalGraphCursor(ReversalGraphCursor&&) = default;
  ReversalGraphCursor& operator=(const ReversalGraphCursor&) = default;
  ReversalGraphCursor& operator=(ReversalGraphCursor&&) = default;

  std::vector<ReversalGraphCursor> next() const {
    auto result = GraphCursor::prev();
    return std::vector<ReversalGraphCursor>(std::cbegin(result), std::cend(result));
  }

  std::vector<ReversalGraphCursor> prev() const {
    auto result = GraphCursor::next();
    return std::vector<ReversalGraphCursor>(std::cbegin(result), std::cend(result));
  }
};

template <class GraphCursor>
class RestrictedGraphCursor : public GraphCursor {
 public:
  RestrictedGraphCursor(const GraphCursor &other, const std::unordered_set<GraphCursor> &space) : GraphCursor(other), space_{space} {}
  RestrictedGraphCursor(GraphCursor &&other, const std::unordered_set<GraphCursor> &space) : GraphCursor(std::move(other)), space_{space} {}

  std::vector<RestrictedGraphCursor> next() const {
    return filter_(GraphCursor::next());
  }

  std::vector<RestrictedGraphCursor> prev() const {
    return filter_(GraphCursor::prev ());
  }

 private:
  const std::unordered_set<GraphCursor> space_;

  using GraphCursor::GraphCursor;
  RestrictedGraphCursor(GraphCursor &&other) : GraphCursor(std::move(other)) {}

  std::vector<RestrictedGraphCursor> filter_(std::vector<GraphCursor> v) {
    std::vector<RestrictedGraphCursor> result;
    for (auto &&cursor : v) {
      if (space_.count(cursor)) {
        result.emplace_back(std::move(cursor));
      }
    }

    return result;
  }
};

namespace std {
template <typename GraphCursor>
struct hash<ReversalGraphCursor<GraphCursor>> : public hash<GraphCursor> {};

template <typename GraphCursor>
struct hash<RestrictedGraphCursor<GraphCursor>> : public hash<GraphCursor> {};
}  // namespace std

// vim: set ts=2 sw=2 et :
