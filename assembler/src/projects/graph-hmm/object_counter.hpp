#pragma once

#include <atomic>

namespace impl {
template <typename T>
void update_max(T &shared, const T &candidate) {
  shared = std::max<T>(shared, candidate);
}

template <typename T>
void update_max(std::atomic<T> &shared, const std::atomic<T> &candidate) {
  shared = std::max<T>(shared, candidate);  // TODO implement a proper atomic max-update
}

template <class T, class Counter>
class ObjectCounter {
 public:
  static size_t object_count_max() { return object_count_max_; }
  static size_t object_count_moved() { return object_count_moved_; }
  static size_t object_count_copied() { return object_count_copied_; }
  static size_t object_count_current() { return object_count_current_; }
  static size_t object_count_constructed() { return object_count_constructed_; }

  ObjectCounter() noexcept { object_count_construct_(); }

  ObjectCounter(const ObjectCounter &) noexcept {
    ++object_count_copied_;
    object_count_construct_();
  }

  ObjectCounter(ObjectCounter &&) noexcept {
    ++object_count_moved_;
    object_count_construct_();
  }

  ~ObjectCounter() { --object_count_current_; }

 private:
  static Counter object_count_max_;
  static Counter object_count_moved_;
  static Counter object_count_copied_;
  static Counter object_count_current_;
  static Counter object_count_constructed_;

  static void object_count_construct_() {
    ++object_count_constructed_;
    ++object_count_current_;
    update_max(object_count_max_, object_count_current_);
  }
};

template <class T, class Counter>
Counter ObjectCounter<T, Counter>::object_count_max_(0);

template <class T, class Counter>
Counter ObjectCounter<T, Counter>::object_count_moved_(0);

template <class T, class Counter>
Counter ObjectCounter<T, Counter>::object_count_copied_(0);

template <class T, class Counter>
Counter ObjectCounter<T, Counter>::object_count_current_(0);

template <class T, class Counter>
Counter ObjectCounter<T, Counter>::object_count_constructed_(0);

}  // namespace impl

template <class T>
using ObjectCounter = impl::ObjectCounter<T, size_t>;

template <class T>
using AtomicObjectCounter = impl::ObjectCounter<T, std::atomic<size_t>>;

class DummyObjectCounter {
 public:
  static size_t object_count_max() { return 0; }
  static size_t object_count_moved() { return 0; }
  static size_t object_count_copied() { return 0; }
  static size_t object_count_current() { return 0; }
  static size_t object_count_constructed() { return 0; }
};

// vim: set ts=2 sw=2 et :
