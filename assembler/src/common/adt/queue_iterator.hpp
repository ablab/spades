//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef QUEUE_ITERATOR_HPP_
#define QUEUE_ITERATOR_HPP_

#include "utils/verify.hpp"
#include "btree/btree_set.h"
#include <set>
#include <queue>
#include <parallel_hashmap/phmap.h>
#include <vector>
#include <algorithm>

namespace {
template<typename T>
inline void hash_combine(std::size_t& seed, const T& val) {
    std::hash<T> hasher;
    seed ^= hasher(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
}

//  taken from https://stackoverflow.com/a/7222201/916549
//
namespace std {
template<typename S, typename T>
struct hash<std::pair<S, T>> {
    inline size_t operator()(const std::pair<S, T>& val) const {
        size_t seed = 0;
        hash_combine(seed, val.first);
        hash_combine(seed, val.second);
        return seed;
    }
};
}

namespace adt {

class identity {
public:
    template <typename T>
    constexpr T&& operator()(T&& t) const noexcept {
        return std::forward<T>(t);
    }
};

// std::greater uses operator> which could be unimplemented
struct simple_greater {
    template <typename T>
    constexpr bool operator()(const T &t1, const T &t2) const noexcept {
        return t2 < t1;
    }
};

template<typename T, typename Priority=identity>
class erasable_priority_queue_key_dirty_heap {
private:
    Priority priority_;
    using PriorityValue = std::decay_t<decltype(std::declval<Priority>()(T()))>;

    constexpr static bool is_trivial = std::is_same<Priority, identity>::value;
    using StoredType = std::conditional_t<is_trivial, T, std::pair<PriorityValue, T>>;

    template <typename... Args>
    class priority_queue : public std::priority_queue<Args...> {
    public:  // use Public Morozov pattern
        using std::priority_queue<Args...>::c;

        // std::priority_queue has not .clear() method!
        void clear() {
            c.clear();
        }
    };

    priority_queue<StoredType, std::vector<StoredType>, simple_greater> queue_;
    phmap::flat_hash_set<StoredType> set_;

    template <typename This>
    struct TrivialFunctions {
        static const T& top_impl(const This &p) {
            return p.queue_.top();
        }

        static const StoredType get_stored(const This&, const T &key) {
            return key;
        }
    };

    template <typename This>
    struct NontrivialFunctions {
        static const T& top_impl(const This &p) {
            return p.queue_.top().second;
        }

        static const StoredType get_stored(const This &p, const T &key) {
            return std::make_pair(p.priority_(key), key);
        }
    };

    using Queue = erasable_priority_queue_key_dirty_heap<T, Priority>;
    using Functions = std::conditional_t<is_trivial, TrivialFunctions<Queue>, NontrivialFunctions<Queue>>;

    const T& top_impl() const {
        return Functions::top_impl(*this);
    }

    StoredType get_stored(const T &key) const {
        return Functions::get_stored(*this, key);
    }

    void skip() {
        if (empty()) return clear();

        while (!set_.count(queue_.top())) {
            queue_.pop();
        }
    }
public:
    erasable_priority_queue_key_dirty_heap() {}
    erasable_priority_queue_key_dirty_heap(Priority priority)
            : priority_(std::move(priority)) {}

    template<typename InputIterator>
    erasable_priority_queue_key_dirty_heap(InputIterator begin, InputIterator end) {
        insert(begin, end);
    }
    template<typename InputIterator>
    erasable_priority_queue_key_dirty_heap(InputIterator begin, InputIterator end,
                                           Priority priority)
            : priority_(std::move(priority)) {
        insert(begin, end);
    }

    void pop() {
        VERIFY(!set_.empty());
        bool res = set_.erase(queue_.top());
        VERIFY(res);
        queue_.pop();
        skip();
    }

    const T& top() const {
        VERIFY(!set_.empty());
        return top_impl();
    }

    void push(const T &key) {
        auto p = get_stored(key);
        queue_.push(p);
        set_.insert(p);
    }

    bool erase(const T &key) {
        auto p = get_stored(key);
        bool res = set_.erase(p) > 0;
        skip();
        if (2 * set_.size() < queue_.size())
            compress();
        return res;
    }

    void clear() {
        set_.clear();
        queue_.clear();
    }

    void compress() {
        auto &storage = queue_.c;
        storage = std::vector<StoredType>(set_.cbegin(), set_.cend());
        std::make_heap(storage.begin(), storage.end(), simple_greater());
    }

    bool empty() const {
        return set_.empty();
    }

    size_t size() const {
        return set_.size();
    }

    template <class InputIterator>
    void insert(InputIterator begin, InputIterator end) {
        for (; begin != end; ++begin) {
            push(*begin);
        }
    }
};

template<typename T, typename Priority=identity>
class erasable_priority_queue_key {
private:
    using PriorityValue = std::decay_t<decltype(std::declval<Priority>()(T()))>;
    btree::btree_set<std::pair<PriorityValue, T>> storage_;
    Priority priority_;
public:
    erasable_priority_queue_key() {}
    erasable_priority_queue_key(Priority priority)
            : priority_(std::move(priority)) {}

    template<typename InputIterator>
    erasable_priority_queue_key(InputIterator begin, InputIterator end) {
        insert(begin, end);
    }
    template<typename InputIterator>
    erasable_priority_queue_key(InputIterator begin, InputIterator end,
                                Priority priority)
            : priority_(std::move(priority)) {
        insert(begin, end);
    }

    void pop() {
        VERIFY(!storage_.empty());
        storage_.erase(storage_.begin());
    }

    const T& top() const {
        VERIFY(!storage_.empty());
        return storage_.begin()->second;
    }

    void push(const T &key) {
        storage_.insert(std::make_pair(priority_(key), key));
    }

    bool erase(const T &key) {
        bool res = storage_.erase(std::make_pair(priority_(key), key)) > 0;
        return res;
    }

    void clear() {
        storage_.clear();
    }

    bool empty() const {
        return storage_.empty();
    }

    size_t size() const {
        return storage_.size();
    }

    template <class InputIterator>
    void insert(InputIterator begin, InputIterator end) {
        for (; begin != end; ++begin) {
            push(*begin);
        }
    }
};

template<typename T, typename Comparator = std::less<T>>
class erasable_priority_queue {
private:
    Comparator cmp_;
    btree::btree_set<T, Comparator> storage_;
public:
    erasable_priority_queue()
            : cmp_(), storage_(cmp_) {}

    erasable_priority_queue(Comparator comparator)
            : cmp_(std::move(comparator)), storage_(cmp_) {}

    template<typename InputIterator>
    erasable_priority_queue(InputIterator begin, InputIterator end)
            : cmp_(), storage_(begin, end, cmp_) {}

    template<typename InputIterator>
    erasable_priority_queue(InputIterator begin, InputIterator end,
                            Comparator comparator)
            : cmp_(std::move(comparator)), storage_(begin, end, cmp_) {}

    void pop() {
        VERIFY(!storage_.empty());
        storage_.erase(storage_.begin());
    }

    const T& top() const {
        VERIFY(!storage_.empty());
        return *(storage_.begin());
    }

    void push(const T &key) {
        storage_.insert(key);
    }

    bool erase(const T &key) {
        bool res = storage_.erase(key) > 0;
        return res;
    }

    void clear() {
        storage_.clear();
    }

    bool empty() const {
        return storage_.empty();
    }

    size_t size() const {
        return storage_.size();
    }

    template <class InputIterator>
    void insert ( InputIterator first, InputIterator last ) {
        storage_.insert(first, last);
    }

};


template<typename T, typename Comparator = std::less<T>>
class indexed_heap_erasable_priority_queue {
private:
    std::vector<T> heap_;
    phmap::flat_hash_map<T, size_t> map_;
    Comparator cmp_;

    template<typename InputIterator>
    void init(InputIterator begin, InputIterator end) {
        for (; begin != end; ++begin) {
            if (map_.insert({*begin, size_t(-1)}).second) {
                heap_.push(*begin);
            }
        }
        auto rev_cmp = [this](const T &a, const T &b) -> bool { return cmp_(b, a); };
        std::make_heap(heap_.begin(), heap_.end(), rev_cmp);
        for (size_t i = 0; i < heap_.size(); ++i) {
            map_[heap_[i]] = i;
        }
        VERIFY(heap_.size() == map_.size());
    }
public:
    indexed_heap_erasable_priority_queue() {}
    indexed_heap_erasable_priority_queue(Comparator comparator = Comparator())
            : cmp_(std::move(comparator)) {}

    template<typename InputIterator>
    indexed_heap_erasable_priority_queue(InputIterator begin, InputIterator end) {
        init(begin, end);
    }

    template<typename InputIterator>
    indexed_heap_erasable_priority_queue(InputIterator begin, InputIterator end,
                                         Comparator comparator)
            : cmp_(std::move(comparator)) {
        init(begin, end);
    }

    void pop() {
        VERIFY_DEV(!heap_.empty());
        heap_pop();
        VERIFY_DEV(map_.size() == heap_.size());
    }

    const T& top() {
        VERIFY_DEV(!heap_.empty());
        return heap_top();
    }

    void push(const T &key) {
        auto it_fl = map_.insert({key, size_t(-1)});
        if (!it_fl.second) return;
        it_fl.first->second = heap_.size();
        heap_.push_back(key);
        up(heap_.size() - 1);
        VERIFY_DEV(map_.size() == heap_.size());
    }

    bool erase(const T &key) {
        auto it = map_.find(key);
        if (it == map_.end()) return false;
        size_t i = it->second;
        heap_[i] = heap_.back();
        map_.find(heap_[i])->second = i;
        map_._erase(it);
        heap_.resize(heap_.size() - 1);
        if (i == heap_.size()) return true;
        if (up(i) == i) down(i);
        VERIFY_DEV(map_.size() == heap_.size());
        return true;
    }

    void clear() {
        map_.clear();
        heap_.clear();
    }

    bool empty() const {
        return heap_.empty();
    }

    size_t size() const {
        return heap_.size();
    }

    template <class InputIterator>
    void insert(InputIterator first, InputIterator last) {
        for (; first != last; ++first) {
            push(*first);
        }
    }

private:
    size_t down(size_t i) {
        T val = heap_[i];
        while (2 * i + 1 < heap_.size()) {
            size_t son = 2 * i + 1;
            size_t other_son = 2 * i + 2;
            if (other_son < heap_.size() && cmp_(heap_[other_son], heap_[son])) son = other_son;
            if (cmp_(heap_[son], val)) {
                heap_[i] = heap_[son];
                map_[heap_[i]] = i;
                i = son;
            } else {
                break;
            }
        }
        heap_[i] = val;
        map_[val] = i;
        return i;
    }

    size_t up(size_t i) {
        T val = heap_[i];
        while (i > 0) {
            size_t parent = (i - 1) / 2;
            if (cmp_(val, heap_[parent])) {
                heap_[i] = heap_[parent];
                map_[heap_[i]] = i;
                i = parent;
            } else {
                break;
            }
        }
        heap_[i] = val;
        map_[val] = i;
        return i;
    }

    void heap_pop() {
        erase(heap_top());
    }

    const T& heap_top() const {
        return heap_[0];
    }
};

template<typename T, typename Queue>
class DynamicQueueIteratorBase {
    bool current_actual_;
    bool current_deleted_;
    T current_;
    Queue queue_;
public:
    template<class... Args>
    DynamicQueueIteratorBase(Args&&... args)
            :  current_actual_(false), current_deleted_(false),
               queue_(std::forward<Args>(args)...) {}

    template<typename InputIterator>
    void insert(InputIterator begin, InputIterator end) {
        queue_.insert(begin, end);
    }

    void push(const T &to_add) {
        queue_.push(to_add);
    }

    void erase(const T &to_remove) {
        if (current_actual_ && to_remove == current_) {
            current_deleted_ = true;
        }
        queue_.erase(to_remove);
    }

    void clear() {
        queue_.clear();
        current_actual_ = false;
        current_deleted_ = false;
    }

    bool IsEnd() const {
        return queue_.empty();
    }

    size_t size() const {
        return queue_.size();
    }

    const T& operator*() {
        VERIFY(!queue_.empty());
        if (!current_actual_ || current_deleted_) {
            current_ = queue_.top();
            current_actual_ = true;
            current_deleted_ = false;
        }
        return current_;
    }

    void operator++() {
        if (!current_actual_) {
            queue_.pop();
        } else if (!current_deleted_) {
            queue_.erase(current_);
        }
        current_actual_ = false;
    }

    //use carefully!
    void ReleaseCurrent() {
        current_actual_ = false;
    }
};

// Iterator over queue that is ordered using Priority
template<typename T, typename Priority = identity>
class DynamicQueueIteratorKey : public DynamicQueueIteratorBase<T, erasable_priority_queue_key_dirty_heap<T, Priority>> {
    using base = DynamicQueueIteratorBase<T, erasable_priority_queue_key_dirty_heap<T, Priority>>;
public:
    DynamicQueueIteratorKey(const Priority &priority = Priority())
            : base(priority) {}
};

// Iterator over quuue that is ordered using Comparator
template<typename T, typename Comparator = std::less<T>>
class DynamicQueueIterator : public DynamicQueueIteratorBase<T, indexed_heap_erasable_priority_queue<T, Comparator>> {
    using base = DynamicQueueIteratorBase<T, indexed_heap_erasable_priority_queue<T, Comparator>>;
public:
    DynamicQueueIterator(const Comparator &comparator = Comparator())
            : base(comparator) {}
};
} //adt

#endif /* QUEUE_ITERATOR_HPP_ */

