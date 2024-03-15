//***************************************************************************
//* Copyright (c) 2018-2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <limits>
#include <cmath>
#include "utils/logger/logger.hpp"
#include <unordered_map>
#include <unordered_set>
#include "utils.hpp"

namespace depth_filter {

namespace impl {

template <typename GraphCursor>
class Depth {
 public:
  bool depth_at_least(const GraphCursor &cursor, double d, typename GraphCursor::Context context) {
    return depth(cursor, context) >= d;
  }

  double depth(const GraphCursor &cursor, typename GraphCursor::Context context) {
    std::unordered_set<GraphCursor> stack;  // TODO do not construct stack in case of using cached value
    assert(stack.size() == 0);
    auto result = get_depth_(cursor, stack, context);
    assert(stack.size() == 0);
    return result;
  }

 size_t max_stack_size() const { return max_stack_size_; }

 private:
  std::unordered_map<GraphCursor, double> depth_;
  size_t max_stack_size_ = 0;

  double get_depth_(const GraphCursor &cursor, std::unordered_set<GraphCursor> &stack, typename GraphCursor::Context context) {
    if (depth_.count(cursor)) {
      return depth_[cursor];
    }

    if (cursor.is_empty()) {
      return depth_[cursor] = std::numeric_limits<double>::infinity();
    }

    if (cursor.letter(context) == '*' || cursor.letter(context) == 'X') {  // FIXME X is not stop codon
      // INFO("Empty depth " << cursor);
      return depth_[cursor] = 0;
    }

    if (stack.count(cursor)) {
      return depth_[cursor] = std::numeric_limits<double>::infinity();
    }

    auto nexts = cursor.next(context);
    stack.insert(cursor);
    max_stack_size_ = std::max(max_stack_size_, stack.size());
    double max_child = 0;
    for (const GraphCursor &n : nexts) {
      max_child = std::max(max_child, get_depth_(n, stack, context));
    }
    stack.erase(cursor);

    return depth_[cursor] = 1 + max_child;
  }
};

template <typename GraphCursor>
class DepthInt {
 public:
  static const size_t INF;

  DepthInt() {
    depth_[GraphCursor()] = INF;
  }

  bool depth_at_least(const GraphCursor &cursor, double d, typename GraphCursor::Context context) {
    return static_cast<double>(depth(cursor, context)) >= d;
  }

  size_t depth(const GraphCursor &cursor, typename GraphCursor::Context context) {
    auto it = depth_.find(cursor);
    if (it != depth_.cend()) {
      return it->second;
    } else {
      return get_depth_non_rec_(cursor, context);
    }
  }


  size_t max_stack_size() const { return max_stack_size_; }

 private:
  std::unordered_map<GraphCursor, size_t> depth_;
  size_t max_stack_size_ = 0;

  // Not thread-safe!
  size_t get_depth_non_rec_(const GraphCursor &cursor, typename GraphCursor::Context context) {
    std::vector<GraphCursor> stack;
    stack.push_back(cursor);

    while (!stack.empty()) {
      max_stack_size_ = std::max(max_stack_size_, stack.size());

      GraphCursor cursor = stack.back();

      if (cursor.letter(context) == '*' || cursor.letter(context) == 'X') {  // FIXME X is not actual stop codon
        depth_[cursor] = 0;
        stack.pop_back();
        continue;
      }

      auto nexts = cursor.next(context);
      // Check children
      size_t max_child = 0;
      size_t unknown_children = 0;
      for (const GraphCursor &n : nexts) {
        auto it = depth_.find(n);
        if (it != depth_.cend()) {
          max_child = std::max(max_child, it->second);
          if (max_child == INF) {
            break;
          }
        } else {
          stack.push_back(n);
          ++unknown_children;
        }
      }

      if (!unknown_children || max_child == INF) {
        stack.resize(stack.size() - unknown_children - 1);
        depth_[cursor] = saturated_inc(max_child, INF);
      } else {
        depth_[cursor] = INF;  // temporary INF to indicate visited vertex
      }
    }

    return depth_.find(cursor)->second;
  }
};

template <typename GraphCursor>
const size_t DepthInt<GraphCursor>::INF = std::numeric_limits<size_t>::max();

template <typename GraphCursor>
class DummyDepthAtLeast {
  bool depth_at_least(const GraphCursor &, double, const void*) {
    return true;
  }
};

template <typename GraphCursor>
class DepthAtLeast {
  struct Estimation {
    size_t value;
    bool exact;
  };

 public:
  static const size_t STACK_LIMIT = 50000;
  static const size_t INF = std::numeric_limits<size_t>::max();

  bool depth_at_least(const GraphCursor &cursor, double depth,
                      typename GraphCursor::Context context) {
    if (depth <= 0) {
      return true;
    }
    return depth_at_least(cursor, static_cast<size_t>(ceil(depth)), context);
  }

  bool depth_at_least(const GraphCursor &cursor, size_t depth,
                      typename GraphCursor::Context context) {
    if (depth == 0) {
      return true;
    }

    if (depth_.count(cursor)) {
      const auto &cached = depth_.find(cursor)->second;
      if (cached.value >= depth) {
        return true;
      } else if (cached.exact) {
        return false;
      }
    }

    const size_t coef = 2;
    size_t stack_limit = std::max<size_t>(coef * depth, 10);

    assert(stack_limit >= depth);
    assert(stack_limit <= STACK_LIMIT);

    std::unordered_set<GraphCursor> stack;
    get_depth_(cursor, stack, stack_limit, context);
    assert(stack.empty());

    assert(depth_.count(cursor));
    return depth_at_least(cursor, depth, context);
  }

  size_t max_stack_size() const { return max_stack_size_; }

 private:
  std::unordered_map<GraphCursor, Estimation> depth_;
  size_t max_stack_size_ = 0;

  Estimation get_depth_(const GraphCursor &cursor, std::unordered_set<GraphCursor> &stack, size_t stack_limit,
                        typename GraphCursor::Context context) {
    if (cursor.is_empty()) {
      return depth_[cursor] = {INF, true};
    }

    if (depth_.count(cursor)) {
      if (depth_[cursor].exact || depth_[cursor].value > stack_limit) {
        return depth_[cursor];
      }
    }

    if (cursor.letter(context) == '*' || cursor.letter(context) == 'X') {  // FIXME X is not stop codon
      return depth_[cursor] = {0, true};
    }

    if (stack.count(cursor)) {
      return depth_[cursor] = {INF, true};
    }

    if (!stack_limit) {
      return depth_[cursor] = {1, false};
    }

    auto nexts = cursor.next(context);
    stack.insert(cursor);
    max_stack_size_ = std::max(max_stack_size_, stack.size());
    size_t max_child = 0;
    bool exact = true;
    for (const GraphCursor &n : nexts) {
      auto result = get_depth_(n, stack, stack_limit - 1, context);
      max_child = std::max(max_child, result.value);
      exact = exact && result.exact;
    }

    stack.erase(cursor);

    if (max_child == INF) {
      // Infinity is always "exact"
      return depth_[cursor] = {INF, true};
    } else {
      return depth_[cursor] = {1 + max_child, exact};
    }
  }
};
} // namespace impl

using impl::DepthInt;
} // namespace depth_filter

// vim: set ts=2 sw=2 et :
