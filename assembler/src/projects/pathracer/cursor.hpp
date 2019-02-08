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

  std::vector<ReversalGraphCursor> next(const void *context) const {
    auto result = GraphCursor::prev(context);
    return std::vector<ReversalGraphCursor>(std::cbegin(result), std::cend(result));
  }

  std::vector<ReversalGraphCursor> prev(const void *context) const {
    auto result = GraphCursor::next(context);
    return std::vector<ReversalGraphCursor>(std::cbegin(result), std::cend(result));
  }
};

template <class GraphCursor>
class RestrictedGraphCursor : public GraphCursor {
 public:
  RestrictedGraphCursor() = default;
  RestrictedGraphCursor(const RestrictedGraphCursor&) = default;
  RestrictedGraphCursor(RestrictedGraphCursor&&) = default;
  RestrictedGraphCursor& operator=(const RestrictedGraphCursor&) = default;
  RestrictedGraphCursor& operator=(RestrictedGraphCursor&&) = default;

  RestrictedGraphCursor(const GraphCursor &other, const std::unordered_set<GraphCursor> &space) : GraphCursor(other), pspace_{&space} {}
  RestrictedGraphCursor(GraphCursor &&other, const std::unordered_set<GraphCursor> &space) : GraphCursor(std::move(other)), pspace_{&space} {}

  std::vector<RestrictedGraphCursor> next(const void *context) const {
    return filter_(GraphCursor::next(context));
  }

  std::vector<RestrictedGraphCursor> prev(const void *context) const {
    return filter_(GraphCursor::prev(context));
  }

 private:
  const std::unordered_set<GraphCursor> *pspace_;

  std::vector<RestrictedGraphCursor> filter_(std::vector<GraphCursor> v) const {
    std::vector<RestrictedGraphCursor> result;
    for (auto &&cursor : v) {
      if (pspace_->count(cursor)) {
        result.emplace_back(std::move(cursor), *pspace_);
      }
    }

    return result;
  }
};

template <class GraphCursor>
auto make_restricted_cursors(const std::vector<GraphCursor> &cursors, const std::unordered_set<GraphCursor> &space) {
  std::vector<RestrictedGraphCursor<GraphCursor>> result;
  for (const auto &cursor : cursors) {
    result.emplace_back(cursor, space);
  }

  return result;
}

template <class GraphCursor>
struct OptimizedRestrictedGraphCursorContext {
    const std::unordered_set<GraphCursor> space;
    const void *context;
};

template <class GraphCursor>
class OptimizedRestrictedGraphCursor : public GraphCursor {
 public:
  OptimizedRestrictedGraphCursor() = default;
  OptimizedRestrictedGraphCursor(const OptimizedRestrictedGraphCursor&) = default;
  OptimizedRestrictedGraphCursor(OptimizedRestrictedGraphCursor&&) = default;
  OptimizedRestrictedGraphCursor& operator=(const OptimizedRestrictedGraphCursor&) = default;
  OptimizedRestrictedGraphCursor& operator=(OptimizedRestrictedGraphCursor&&) = default;

  OptimizedRestrictedGraphCursor(const GraphCursor &other) : GraphCursor(other) {}
  OptimizedRestrictedGraphCursor(GraphCursor &&other) : GraphCursor(std::move(other)) {}

  std::vector<OptimizedRestrictedGraphCursor> next(const void *context) const {
    return filter_(GraphCursor::next(unvoid_context(context)->context), unvoid_context(context)->space);
  }

  std::vector<OptimizedRestrictedGraphCursor> prev(const void *context) const {
    return filter_(GraphCursor::prev(unvoid_context(context)->context), unvoid_context(context)->space);
  }

  char letter(const void *context) const { return GraphCursor::letter(unvoid_context(context)->context); }

 private:
  static const OptimizedRestrictedGraphCursorContext<GraphCursor> *unvoid_context(const void *context) {
    return static_cast<const OptimizedRestrictedGraphCursorContext<GraphCursor>*>(context);
  }

  std::vector<OptimizedRestrictedGraphCursor> filter_(std::vector<GraphCursor> v, const std::unordered_set<GraphCursor> &space) const {
    std::vector<OptimizedRestrictedGraphCursor> result;
    for (auto &&cursor : v) {
      if (space.count(cursor)) {
        result.emplace_back(std::move(cursor));
      }
    }

    return result;
  }
};

template <class GraphCursor>
auto make_optimized_restricted_cursor_context(const std::unordered_set<GraphCursor> &space, const void *context) {
    return OptimizedRestrictedGraphCursorContext<GraphCursor>{space, context};
}

template <class GraphCursor>
auto make_optimized_restricted_cursors(const std::vector<GraphCursor> &cursors) {
  std::vector<OptimizedRestrictedGraphCursor<GraphCursor>> result;
  for (const auto &cursor : cursors) {
    result.emplace_back(cursor);
  }

  return result;

}
namespace std {
template <typename GraphCursor>
struct hash<ReversalGraphCursor<GraphCursor>> : public hash<GraphCursor> {};

template <typename GraphCursor>
struct hash<RestrictedGraphCursor<GraphCursor>> : public hash<GraphCursor> {};

template <typename GraphCursor>
struct hash<OptimizedRestrictedGraphCursor<GraphCursor>> : public hash<GraphCursor> {};
}  // namespace std

// vim: set ts=2 sw=2 et :
