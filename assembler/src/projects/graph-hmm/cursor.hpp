#pragma once

#include <vector>
#include <utility>

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

namespace std {
template <typename GraphCursor>
struct hash<ReversalGraphCursor<GraphCursor>> : public hash<GraphCursor> {};
}  // namespace std
