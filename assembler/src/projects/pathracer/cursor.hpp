#pragma once

#include <vector>
#include <utility>
#include <unordered_set>

template <class GraphCursor>
struct OptimizedRestrictedGraphCursorContext {
    const std::unordered_set<GraphCursor> &space;
    typename GraphCursor::Context context;
};

template <class GraphCursor>
class OptimizedRestrictedGraphCursor : public GraphCursor {
 public:
  using Context = const OptimizedRestrictedGraphCursorContext<GraphCursor>*;
  OptimizedRestrictedGraphCursor() = default;
  OptimizedRestrictedGraphCursor(const OptimizedRestrictedGraphCursor&) = default;
  OptimizedRestrictedGraphCursor(OptimizedRestrictedGraphCursor&&) = default;
  OptimizedRestrictedGraphCursor& operator=(const OptimizedRestrictedGraphCursor&) = default;
  OptimizedRestrictedGraphCursor& operator=(OptimizedRestrictedGraphCursor&&) = default;

  OptimizedRestrictedGraphCursor(const GraphCursor &other) : GraphCursor(other) {}
  OptimizedRestrictedGraphCursor(GraphCursor &&other) : GraphCursor(std::move(other)) {}

  std::vector<OptimizedRestrictedGraphCursor> next(Context context) const {
    return filter_(GraphCursor::next(context->context), context->space);
  }

  std::vector<OptimizedRestrictedGraphCursor> prev(Context context) const {
    return filter_(GraphCursor::prev(context->context), context->space);
  }

  char letter(Context context) const { return GraphCursor::letter(context->context); }

 private:
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
auto make_optimized_restricted_cursor_context(const std::unordered_set<GraphCursor> &space, typename GraphCursor::Context context) {
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
struct hash<OptimizedRestrictedGraphCursor<GraphCursor>> : public hash<GraphCursor> {};
}  // namespace std

// vim: set ts=2 sw=2 et :
