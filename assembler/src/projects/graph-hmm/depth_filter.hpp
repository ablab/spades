#pragma once

#include <limits>
#include <cmath>
#include "utils/logger/logger.hpp"
#include <unordered_map>
#include <unordered_set>

namespace impl {

template <typename GraphCursor>
class Depth {
 public:
  bool depth_at_least(const GraphCursor &cursor, double d) {
    return depth(cursor) >= d;
  }

  double depth(const GraphCursor &cursor) {
    std::unordered_set<GraphCursor> stack;  // TODO do not construct stack in case of using cached value
    assert(stack.size() == 0);
    auto result = get_depth_(cursor, stack);
    assert(stack.size() == 0);
    return result;
  }

 size_t max_stack_size() const { return max_stack_size_; }

 private:
  std::unordered_map<GraphCursor, double> depth_;
  size_t max_stack_size_ = 0;

  double get_depth_(const GraphCursor &cursor, std::unordered_set<GraphCursor> &stack) {
    if (depth_.count(cursor)) {
      return depth_[cursor];
    }

    if (cursor.is_empty()) {
      return depth_[cursor] = std::numeric_limits<double>::infinity();
    }

    if (cursor.letter() == '*') {
      // INFO("Empty depth " << cursor);
      return depth_[cursor] = 0;
    }

    if (stack.count(cursor)) {
      return depth_[cursor] = std::numeric_limits<double>::infinity();
    }

    auto nexts = cursor.next();
    stack.insert(cursor);
    max_stack_size_ = std::max(max_stack_size_, stack.size());
    double max_child = 0;
    for (const GraphCursor &n : nexts) {
      max_child = std::max(max_child, get_depth_(n, stack));
    }
    stack.erase(cursor);

    return depth_[cursor] = 1 + max_child;
  }
};

template <typename GraphCursor>
class DepthAtLeast {
  struct Estimation {
    double value;
    bool exact;
  };

 public:
  bool depth_at_least(const GraphCursor &cursor, double d) {
    if (d <= 0) {
      return true;
    }

    if (depth_.count(cursor)) {
      const auto &cached = depth_.find(cursor)->second;
      if (cached.value >= d) {
        return true;
      } else if (cached.exact) {
        return false;
      }
    }

    std::unordered_set<GraphCursor> stack;
    const double coef = 2.0;
    size_t stack_limit = static_cast<size_t>(coef * d);
    stack_limit = std::max<size_t>(stack_limit, 10);
    get_depth_(cursor, stack, stack_limit);

    assert(depth_.count(cursor));
    return depth_at_least(cursor, d);
  }

  size_t max_stack_size() const { return max_stack_size_; }

 private:
  std::unordered_map<GraphCursor, Estimation> depth_;
  size_t max_stack_size_ = 0;

  Estimation get_depth_(const GraphCursor &cursor, std::unordered_set<GraphCursor> &stack, size_t stack_limit) {
    assert(stack_limit < 100000);

    if (cursor.is_empty()) {
      return depth_[cursor] = {std::numeric_limits<double>::infinity(), true};
    }

    if (cursor.letter() == '*') {
      return depth_[cursor] = {0, true};
    }

    if (stack.count(cursor)) {
      return depth_[cursor] = {std::numeric_limits<double>::infinity(), true};
    }

    if (stack.size() > stack_limit) {
      return depth_[cursor] = {1, false};
    }

    auto nexts = cursor.next();
    stack.insert(cursor);
    max_stack_size_ = std::max(max_stack_size_, stack.size());
    double max_child = 0;
    bool exact = true;
    for (const GraphCursor &n : nexts) {
      auto result = get_depth_(n, stack, stack_limit);
      max_child = std::max(max_child, result.value);
      exact = exact && result.exact;
    }

    // Infinity is always "exact"
    if (std::isinf(max_child)) {
      exact = true;
    }

    stack.erase(cursor);

    return depth_[cursor] = {1 + max_child, exact};
  }
};

}  // namespace impl
