#pragma once

#include <vector>
#include <utility>

class EmptyClass {};

template <class GraphCursor, class Base = EmptyClass>
class AbstractGraphCursor : public Base {
 public:
  AbstractGraphCursor() = default;
  AbstractGraphCursor(const Base &gc) : Base{gc} {}
  AbstractGraphCursor(Base &&gc) : Base{gc} {}

  bool is_divergent() const { return crtp_this()->next().size() > 1; }
  bool is_convergent() const { return crtp_this()->prev().size() > 1; }

  std::vector<std::pair<GraphCursor, char>> next_pairs() const {
    std::vector<std::pair<GraphCursor, char>> result;
    for (const auto cur : crtp_this()->next()) {
      result.emplace_back(cur, cur.letter());
    }
    return result;
  }

  std::vector<std::pair<GraphCursor, char>> prev_pairs() const {
    std::vector<std::pair<GraphCursor, char>> result;
    for (const auto cur : crtp_this()->prev()) {
      result.emplace_back(cur, cur.letter());
    }
    return result;
  }

 private:
  GraphCursor *crtp_this() { return static_cast<GraphCursor *>(this); }
  const GraphCursor *crtp_this() const { return static_cast<const GraphCursor *>(this); }
};

template <class GraphCursor>
class ReversalGraphCursor : public AbstractGraphCursor<ReversalGraphCursor<GraphCursor>, GraphCursor> {
 public:
  using AbstractGraphCursor<ReversalGraphCursor<GraphCursor>, GraphCursor>::AbstractGraphCursor;
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
