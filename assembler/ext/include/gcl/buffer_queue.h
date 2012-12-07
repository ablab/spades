// Copyright 2010 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef BUFFER_QUEUE_H
#define BUFFER_QUEUE_H

#include <mutex>
#include <condition_variable>

#include "queue_base.h"

namespace gcl {

template <typename Value>
class buffer_queue
{
  public:
    typedef Value value_type;

    buffer_queue() CXX0X_DELETED
    buffer_queue(const buffer_queue&) CXX0X_DELETED
    buffer_queue(size_t max_elems, const char* name);
    explicit buffer_queue(size_t max_elems);
    template <typename Iter>
    buffer_queue(size_t max_elems, Iter first, Iter last, const char* name);
    template <typename Iter>
    buffer_queue(size_t max_elems, Iter first, Iter last);
    buffer_queue& operator =(const buffer_queue&) CXX0X_DELETED
    ~buffer_queue();

//TODO(crowl): Do we want this?
#if 0
    generic_queue_front<value_type> front()
        { return generic_queue_front<value_type>(this); }
    generic_queue_back<value_type> back()
        { return generic_queue_back<value_type>(this); }
#endif

    void close();
    bool is_closed();
    bool is_empty();

    Value value_pop();
    queue_op_status try_pop(Value&);
    queue_op_status wait_pop(Value&);

    void push(const Value& x);
    queue_op_status try_push(const Value& x);
    queue_op_status wait_push(const Value& x);
#ifdef HAS_CXX0X_RVREF
    void push(Value&& x);
    queue_op_status try_push(Value&& x);
    queue_op_status wait_push(Value&& x);
#endif

    const char* name();

  private:
    mutex mtx_;
    condition_variable not_empty_;
    condition_variable not_full_;
    size_t waiting_full_;
    size_t waiting_empty_;
    Value* buffer_;
    size_t push_index_;
    size_t pop_index_;
    size_t num_slots_;
    bool closed_;
    const char* name_;

    void init(size_t max_elems);

    template <typename Iter>
    void iter_init(size_t max_elems, Iter first, Iter last);

    size_t next(size_t idx) { return (idx + 1) % num_slots_; }

    void pop_from(Value& elem, size_t pdx, size_t hdx)
    {
        pop_index_ = next( pdx );
#ifdef HAS_CXX0X_RVREF
        elem = std::move(buffer_[pdx]);
#else
        elem = buffer_[pdx];
#endif
        if ( waiting_full_ > 0 ) {
            --waiting_full_;
            not_full_.notify_one();
        }
    }

    void push_at(const Value& elem, size_t hdx, size_t nxt, size_t pdx)
    {
        buffer_[hdx] = elem;
        push_index_ = nxt;
        if ( waiting_empty_ > 0 ) {
            --waiting_empty_;
            not_empty_.notify_one();
        }
    }

#ifdef HAS_CXX0X_RVREF
    void push_at(Value&& elem, size_t hdx, size_t nxt, size_t pdx)
    {
        buffer_[hdx] = std::move(elem);
        push_index_ = nxt;
        if ( waiting_empty_ > 0 ) {
            --waiting_empty_;
            not_empty_.notify_one();
        }
    }
#endif

};

template <typename Value>
void buffer_queue<Value>::init(size_t max_elems)
{
    if ( max_elems < 1 )
        throw std::invalid_argument("number of elements must be at least one");
}

template <typename Value>
buffer_queue<Value>::buffer_queue(size_t max_elems, const char* name)
:
    waiting_full_( 0 ),
    waiting_empty_( 0 ),
    buffer_( new Value[max_elems+1] ),
    push_index_( 0 ),
    pop_index_( 0 ),
    num_slots_( max_elems+1 ),
    closed_( false ),
    name_( name )
{
    init(max_elems);
}

template <typename Value>
buffer_queue<Value>::buffer_queue(size_t max_elems)
:
    // would rather do buffer_queue(max_elems, "")
    waiting_full_( 0 ),
    waiting_empty_( 0 ),
    buffer_( new Value[max_elems+1] ),
    push_index_( 0 ),
    pop_index_( 0 ),
    num_slots_( max_elems+1 ),
    closed_( false ),
    name_( "" )
{
    init(max_elems);
}

template <typename Value>
template <typename Iter>
void buffer_queue<Value>::iter_init(size_t max_elems, Iter first, Iter last)
{
    size_t hdx = 0;
    for ( Iter cur = first; cur != last; ++cur ) {
        if ( hdx >= max_elems )
            throw std::invalid_argument("too few slots for iterator");
        size_t nxt = hdx + 1; // more efficient than next(hdx)
        size_t pdx = pop_index_;
        push_at( *cur, hdx, nxt, pdx );
        hdx = nxt;
    }
}

template <typename Value>
template <typename Iter>
buffer_queue<Value>::buffer_queue(size_t max_elems, Iter first, Iter last,
                                    const char* name)
:
    // would rather do buffer_queue(max_elems, name)
    waiting_full_( 0 ),
    waiting_empty_( 0 ),
    buffer_( new Value[max_elems+1] ),
    push_index_( 0 ),
    pop_index_( 0 ),
    num_slots_( max_elems+1 ),
    closed_( false ),
    name_( name )
{
    iter_init(max_elems, first, last);
}

template <typename Value>
template <typename Iter>
buffer_queue<Value>::buffer_queue(size_t max_elems, Iter first, Iter last)
:
    // would rather do buffer_queue(max_elems, first, last, "")
    waiting_full_( 0 ),
    waiting_empty_( 0 ),
    buffer_( new Value[max_elems+1] ),
    push_index_( 0 ),
    pop_index_( 0 ),
    num_slots_( max_elems+1 ),
    closed_( false ),
    name_( "" )
{
    iter_init(max_elems, first, last);
}

template <typename Value>
buffer_queue<Value>::~buffer_queue()
{
    delete[] buffer_;
}

template <typename Value>
void buffer_queue<Value>::close()
{
    lock_guard<mutex> hold( mtx_ );
    closed_ = true;
    not_empty_.notify_all();
    not_full_.notify_all();
}

template <typename Value>
bool buffer_queue<Value>::is_closed()
{
    lock_guard<mutex> hold( mtx_ );
    return closed_;
}

template <typename Value>
bool buffer_queue<Value>::is_empty()
{
    lock_guard<mutex> hold( mtx_ );
    return push_index_ == pop_index_;
}

template <typename Value>
queue_op_status buffer_queue<Value>::try_pop(Value& elem)
{
    /* This try block is here to catch exceptions from the mutex
       operations or from the user-defined copy assignment operator
       in the pop_from operation. */
    try {
        lock_guard<mutex> hold( mtx_ );
        size_t pdx = pop_index_;
        size_t hdx = push_index_;
        if ( pdx == hdx ) {
            if ( closed_ )
                return CXX0X_ENUM_QUAL(queue_op_status)closed;
            else
                return CXX0X_ENUM_QUAL(queue_op_status)empty;
        }
        pop_from( elem, pdx, hdx );
        return CXX0X_ENUM_QUAL(queue_op_status)success;
    } catch (...) {
        close();
        throw;
    }
}

template <typename Value>
queue_op_status buffer_queue<Value>::wait_pop(Value& elem)
{
    /* This try block is here to catch exceptions from the mutex
       operations or from the user-defined copy assignment operator
       in the pop_from operation. */
    try {
        unique_lock<mutex> hold( mtx_ );
        size_t pdx;
        size_t hdx;
        for (;;) {
            pdx = pop_index_;
            hdx = push_index_;
            if ( pdx != hdx )
                break;
            if ( closed_ )
                return CXX0X_ENUM_QUAL(queue_op_status)closed;
            ++waiting_empty_;
            not_empty_.wait( hold );
        }
        pop_from( elem, pdx, hdx );
        return CXX0X_ENUM_QUAL(queue_op_status)success;
    } catch (...) {
        close();
        throw;
    }
}

template <typename Value>
Value buffer_queue<Value>::value_pop()
{
    /* This try block is here to catch exceptions from the
       user-defined copy assignment operator. */
    try {
        Value elem;
        if ( wait_pop( elem ) == CXX0X_ENUM_QUAL(queue_op_status)closed )
            throw CXX0X_ENUM_QUAL(queue_op_status)closed;
#ifdef HAS_CXX0X_RVREF
        return std::move(elem);
#else
        return elem;
#endif
    } catch (...) {
        close();
        throw;
    }
}

template <typename Value>
queue_op_status buffer_queue<Value>::try_push(const Value& elem)
{
    /* This try block is here to catch exceptions from the mutex
       operations or from the user-defined copy assignment
       operator in push_at. */
    try {
        lock_guard<mutex> hold( mtx_ );
        if ( closed_ )
            return CXX0X_ENUM_QUAL(queue_op_status)closed;
        size_t hdx = push_index_;
        size_t nxt = next( hdx );
        size_t pdx = pop_index_;
        if ( nxt == pdx )
            return CXX0X_ENUM_QUAL(queue_op_status)full;
        push_at( elem, hdx, nxt, pdx );
        return CXX0X_ENUM_QUAL(queue_op_status)success;
    } catch (...) {
        close();
        throw;
    }
}

template <typename Value>
queue_op_status buffer_queue<Value>::wait_push(const Value& elem)
{
    /* This try block is here to catch exceptions from the mutex
       operations or from the user-defined copy assignment
       operator in push_at. */
    try {
        unique_lock<mutex> hold( mtx_ );
        size_t hdx;
        size_t nxt;
        size_t pdx;
        for (;;) {
            if ( closed_ )
                return CXX0X_ENUM_QUAL(queue_op_status)closed;
            hdx = push_index_;
            nxt = next( hdx );
            pdx = pop_index_;
            if ( nxt != pdx )
                break;
            ++waiting_full_;
            not_full_.wait( hold );
        }
        push_at( elem, hdx, nxt, pdx );
        return CXX0X_ENUM_QUAL(queue_op_status)success;
    } catch (...) {
        close();
        throw;
    }
}

template <typename Value>
void buffer_queue<Value>::push(const Value& elem)
{
    /* Only wait_push can throw, and it protects itself, so there
       is no need to try/catch here. */
    if ( wait_push( elem ) == CXX0X_ENUM_QUAL(queue_op_status)closed ) {
        throw CXX0X_ENUM_QUAL(queue_op_status)closed;
    }
}

#ifdef HAS_CXX0X_RVREF

//TODO(crowl) Refactor with non-move versions.

template <typename Value>
queue_op_status buffer_queue<Value>::try_push(Value&& elem)
{
    /* This try block is here to catch exceptions from the mutex
       operations or from the user-defined copy assignment
       operator in push_at. */
    try {
        lock_guard<mutex> hold( mtx_ );
        if ( closed_ )
            return CXX0X_ENUM_QUAL(queue_op_status)closed;
        size_t hdx = push_index_;
        size_t nxt = next( hdx );
        size_t pdx = pop_index_;
        if ( nxt == pdx )
            return CXX0X_ENUM_QUAL(queue_op_status)full;
        push_at( std::move(elem), hdx, nxt, pdx );
        return CXX0X_ENUM_QUAL(queue_op_status)success;
    } catch (...) {
        close();
        throw;
    }
}

template <typename Value>
queue_op_status buffer_queue<Value>::wait_push(Value&& elem)
{
    /* This try block is here to catch exceptions from the mutex
       operations or from the user-defined copy assignment
       operator in push_at. */
    try {
        unique_lock<mutex> hold( mtx_ );
        size_t hdx;
        size_t nxt;
        size_t pdx;
        for (;;) {
            if ( closed_ )
                return CXX0X_ENUM_QUAL(queue_op_status)closed;
            hdx = push_index_;
            nxt = next( hdx );
            pdx = pop_index_;
            if ( nxt != pdx )
                break;
            ++waiting_full_;
            not_full_.wait( hold );
        }
        push_at( std::move(elem), hdx, nxt, pdx );
        return CXX0X_ENUM_QUAL(queue_op_status)success;
    } catch (...) {
        close();
        throw;
    }
}

template <typename Value>
void buffer_queue<Value>::push(Value&& elem)
{
    /* Only wait_push can throw, and it protects itself, so there
       is no need to try/catch here. */
    if ( wait_push( std::move(elem) )
         == CXX0X_ENUM_QUAL(queue_op_status)closed )
        throw CXX0X_ENUM_QUAL(queue_op_status)closed;
}

#endif

template <typename Value>
const char* buffer_queue<Value>::name()
{
    return name_;
}

} // namespace gcl

#endif
