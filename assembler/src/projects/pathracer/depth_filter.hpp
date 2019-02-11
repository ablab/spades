#pragma once

#include <limits>
#include <cmath>
#include "utils/logger/logger.hpp"
#include <boost/heap/binomial_heap.hpp>
#include <unordered_map>
#include <unordered_set>

namespace depth_filter {

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

    if (cursor.letter() == '*' || cursor.letter() == 'X') {
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
class DepthInt {
 public:
  bool depth_at_least(const GraphCursor &cursor, size_t d) {
    return depth(cursor) >= d;
  }

  size_t depth(const GraphCursor &cursor) {
    std::unordered_set<GraphCursor> stack;  // TODO do not construct stack in case of using cached value
    assert(stack.size() == 0);
    auto result = get_depth_(cursor, stack);
    assert(stack.size() == 0);
    return result;
  }

  size_t max_stack_size() const { return max_stack_size_; }
  static const size_t INF = std::numeric_limits<size_t>::max();

 private:
  std::unordered_map<GraphCursor, size_t> depth_;
  size_t max_stack_size_ = 0;

  size_t get_depth_(const GraphCursor &cursor, std::unordered_set<GraphCursor> &stack) {
    if (depth_.count(cursor)) {
      return depth_[cursor];
    }

    if (cursor.is_empty()) {
      return depth_[cursor] = INF;
    }

    if (cursor.letter() == '*' || cursor.letter() == 'X') {
      // INFO("Empty depth " << cursor);
      return depth_[cursor] = 0;
    }

    if (stack.count(cursor)) {
      return depth_[cursor] = INF;
    }

    auto nexts = cursor.next();
    stack.insert(cursor);
    max_stack_size_ = std::max(max_stack_size_, stack.size());
    size_t max_child = 0;
    for (const GraphCursor &n : nexts) {
      max_child = std::max(max_child, get_depth_(n, stack));
    }
    stack.erase(cursor);

    return depth_[cursor] = (max_child == INF) ? INF : 1 + max_child;
  }
};


template <typename GraphCursor>
class DummyDepthAtLeast {
  bool depth_at_least(const GraphCursor &, double) {
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
                      const void *context) {
    if (depth <= 0) {
      return true;
    }
    return depth_at_least(cursor, static_cast<size_t>(depth), context);
  }

  bool depth_at_least(const GraphCursor &cursor, size_t depth,
                      const void *context) {
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
                        const void *context) {
    if (cursor.is_empty()) {
      return depth_[cursor] = {INF, true};
    }

    if (depth_.count(cursor)) {
      if (depth_[cursor].exact || depth_[cursor].value > stack_limit) {
        return depth_[cursor];
      }
    }

    if (cursor.letter(context) == '*' || cursor.letter(context) == 'X') {
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

template <class Key, class Value>
class IndexedPriorityQueue {
    struct KeyValue {
        Key key;
        Value value;
        bool operator<(const KeyValue &other) const { return value < other.value; }
    };

    using PriorityQueue = boost::heap::binomial_heap<KeyValue>;

public:
    void push_or_increase(const Key &k, const Value &v) {
        if (has(k)) {
            increase(k, v);
        } else {
            push(k, v);
        }
    }

    void push(const Key &k, const Value &v) { queue_.push({k, v}); }

    void increase(const Key &k, const Value &v) {
        auto handle = index_[k];
        KeyValue kv{k, v};
        queue_.increase(handle, kv);
    }

    bool has(const Key &k) const { return index_.count(k); }

    auto top() const { return queue_.top(); }

    void pop() { queue_.pop(); }

    bool empty() const { return queue_.empty(); }

    size_t size() const { return queue_.size(); }

private:
    PriorityQueue queue_;
    std::unordered_map<Key, typename PriorityQueue::handle_type> index_;
};

template <typename GraphCursor>
std::vector<GraphCursor> depth_subset(const std::vector<std::pair<GraphCursor, size_t>> &initial,
                                      const void *context,
                                      bool forward = true) {
  IndexedPriorityQueue<GraphCursor, size_t> q;

  std::unordered_map<GraphCursor, size_t> initial_map;
  for (const auto &cursor_with_depth : initial) {
    initial_map[cursor_with_depth.first] = std::max(initial_map[cursor_with_depth.first], cursor_with_depth.second);
  }

  for (auto &&cursor_with_depth : initial_map) {
    q.push(cursor_with_depth.first, cursor_with_depth.second);
    // q.push(std::move(cursor_with_depth.first), std::move(cursor_with_depth.second));
  }
  initial_map.clear();

  INFO("Initial queue size: " << q.size());
  std::unordered_set<GraphCursor> visited;
  while (!q.empty()) {
    auto cursor_with_depth = q.top();
    q.pop();

    if (visited.count(cursor_with_depth.key)) {
      continue;
    }
    visited.insert(cursor_with_depth.key);

    if (cursor_with_depth.value > 0) {
      auto cursors = forward ? cursor_with_depth.key.next(context) : cursor_with_depth.key.prev(context);
      for (const auto &cursor : cursors) {
        if (!visited.count(cursor)) {
          q.push_or_increase(cursor, cursor_with_depth.value - 1);
        }
      }
    }
  }

  return std::vector<GraphCursor>(visited.cbegin(), visited.cend());
}

}  // namespace impl

template <typename GraphCursor>
std::vector<GraphCursor> subset(const GraphCursor& cursor, size_t depth,
                                const void *context, bool forward = true) {
    return impl::depth_subset(std::vector<std::pair<GraphCursor, size_t>>({std::make_pair(cursor, depth)}), context, forward);
}

template <typename GraphCursor>
std::vector<GraphCursor> subset(const std::vector<std::pair<GraphCursor, size_t>> &initial,
                                const void *context,
                                bool forward = true) {
    return impl::depth_subset(initial, context, forward);
}

} // namespace depth_filter 

// vim: set ts=2 sw=2 et :
