#pragma once

#include <limits>
#include "utils/logger/logger.hpp"
#include <unordered_map>
#include <unordered_set>

namespace impl {

template <typename GraphCursor>
class Depth {
 public:
  double depth(const GraphCursor &cursor) {
    assert(stack_.size() == 0);
    auto result = get_depth_(cursor);
    assert(stack_.size() == 0);
    return result;
  }

 size_t max_stack_size() const { return max_stack_size_; }

 private:
  std::unordered_map<GraphCursor, double> depth_;
  std::unordered_set<GraphCursor> stack_;
  size_t max_stack_size_ = 0;

  double get_depth_(const GraphCursor &cursor) {
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

    if (stack_.count(cursor)) {
      return depth_[cursor] = std::numeric_limits<double>::infinity();
    }

    auto nexts = cursor.next();
    stack_.insert(cursor);
    max_stack_size_ = std::max(max_stack_size_, stack_.size());
    double max_child = 0;
    for (const GraphCursor &n : nexts) {
      max_child = std::max(max_child, get_depth_(n));
    }
    stack_.erase(cursor);

    return depth_[cursor] = 1 + max_child;
  }
};

}  // namespace impl
