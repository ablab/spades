//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef QUEUE_ITERATOR_HPP_
#define QUEUE_ITERATOR_HPP_

#include "utils/verify.hpp"
#include <set>

template<typename T, typename Comparator>
class erasable_priority_queue {
private:
    std::set<T, Comparator> storage_;
public:
    /*
     * Be careful! This constructor requires Comparator to have default constructor even if you call it with
     * specified comparator. In this case just create default constructor with VERIFY(false) inside it.
     */
    erasable_priority_queue(const Comparator& comparator = Comparator()) :
        storage_(comparator) {
    }

    template<typename InputIterator>
    erasable_priority_queue(InputIterator begin, InputIterator end,
            const Comparator& comparator = Comparator()) :
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

    void push(const T& key) {
        storage_.insert(key);
    }

    bool erase(const T& key) {
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
class DynamicQueueIterator {

    bool current_actual_;
    bool current_deleted_;
    T current_;
    erasable_priority_queue<T, Comparator> queue_;

public:

    DynamicQueueIterator(const Comparator& comparator = Comparator()) :
        current_actual_(false), current_deleted_(false), queue_(comparator) {
    }

    template<typename InputIterator>
    void insert(InputIterator begin, InputIterator end) {
        queue_.insert(begin, end);
    }

    void push(const T& to_add) {
        queue_.push(to_add);
    }

    void erase(const T& to_remove) {
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
        if(!current_actual_ || current_deleted_) {
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


#endif /* QUEUE_ITERATOR_HPP_ */

