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

namespace adt {

class identity {
public:
    template <typename T>
    constexpr T&& operator()(T&& t) const noexcept {
        return std::forward<T>(t);
    }
};

template<typename T, typename Priority=identity>
class erasable_priority_queue_key_dirty_heap {
private:
    using PriorityValue = std::decay_t<decltype(std::declval<Priority>()(T()))>;
    using Pair = std::pair<PriorityValue, T>;
    std::priority_queue<Pair, std::vector<Pair>, std::greater<Pair>> queue_;
    phmap::flat_hash_set<T> set_;
    Priority priority_;

    const T& top_impl() const {
        return queue_.top().second;
    }

    void skip() {
        while (!queue_.empty() && !set_.count(top_impl())) {
            queue_.pop();
        }
    }
public:
    /*
     * Be careful! This constructor requires Priority to have default constructor even if you call it with
     * specified priority. In this case just create default constructor with VERIFY(false) inside it.
     * TODO Fix it
     */
    erasable_priority_queue_key_dirty_heap(const Priority &priority = Priority()) : priority_(priority) {}

    template<typename InputIterator>
    erasable_priority_queue_key_dirty_heap(InputIterator begin, InputIterator end,
                                           const Priority &priority) : priority_(priority) {
        insert(begin, end);
    }

    void pop() {
        VERIFY(!set_.empty());
        set_.erase(top());
        queue_.pop();
        skip();
    }

    const T& top() const {
        VERIFY(!set_.empty());
        return top_impl();
    }

    void push(const T &key) {
        queue_.push(std::make_pair(priority_(key), key));
        set_.insert(key);
    }

    bool erase(const T &key) {
        bool res = set_.erase(key) > 0;
        skip();
        return res;
    }

    void clear() {
        set_.clear();
        queue_ = std::decay_t<decltype(queue_)>();  // std::priority_queue has not .clear() method!
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
    /*
     * Be careful! This constructor requires Priority to have default constructor even if you call it with
     * specified priority. In this case just create default constructor with VERIFY(false) inside it.
     * TODO Fix it
     */
    erasable_priority_queue_key(const Priority &priority = Priority()) : priority_(priority) {}

    template<typename InputIterator>
    erasable_priority_queue_key(InputIterator begin, InputIterator end,
                                const Priority &priority) : priority_(priority) {
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

template<typename T, typename Comparator>
class erasable_priority_queue {
private:
    btree::btree_set<T, Comparator> storage_;
public:
    /*
     * Be careful! This constructor requires Comparator to have default constructor even if you call it with
     * specified comparator. In this case just create default constructor with VERIFY(false) inside it.
     */
    erasable_priority_queue(const Comparator &comparator = Comparator()) :
        storage_(comparator) {
    }

    template<typename InputIterator>
    erasable_priority_queue(InputIterator begin, InputIterator end,
            const Comparator &comparator = Comparator()) :
        storage_(begin, end, comparator) {
    }

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


template<typename T, typename Comparator>
class indexed_heap_erasable_priority_queue {
private:
    std::vector<T> heap_;
    phmap::flat_hash_map<T, size_t> map_;
    Comparator cmp_;
public:
    /*
     * Be careful! This constructor requires Comparator to have default constructor even if you call it with
     * specified comparator. In this case just create default constructor with VERIFY(false) inside it.
     */
    indexed_heap_erasable_priority_queue(const Comparator &comparator = Comparator()) :
        cmp_(comparator) {
    }

    template<typename InputIterator>
    indexed_heap_erasable_priority_queue(InputIterator begin, InputIterator end,
            const Comparator &comparator = Comparator()) :
        cmp_(comparator) {
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

template<typename T, typename Priority = identity>
class DynamicQueueIteratorKey {

    bool current_actual_;
    bool current_deleted_;
    T current_;
    erasable_priority_queue_key_dirty_heap<T, Priority> queue_;

public:

    DynamicQueueIteratorKey(const Priority &priority = Priority()) :
        current_actual_(false), current_deleted_(false), queue_(priority) {
    }

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

template<typename T, typename Comparator = std::less<T>>
class DynamicQueueIterator {

    bool current_actual_;
    bool current_deleted_;
    T current_;
    indexed_heap_erasable_priority_queue<T, Comparator> queue_;

public:

    DynamicQueueIterator(const Comparator &comparator = Comparator()) :
        current_actual_(false), current_deleted_(false), queue_(comparator) {
    }

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
} //adt

#endif /* QUEUE_ITERATOR_HPP_ */

