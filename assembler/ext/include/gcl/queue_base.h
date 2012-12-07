// Copyright 2011 Google Inc. All Rights Reserved.
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

#ifndef QUEUE_BASE_H
#define QUEUE_BASE_H

#include <stddef.h>

#include <iterator>
#include <iostream>

#include "cxx0x.h"
#include "atomic.h"

namespace gcl {

template <typename Queue>
class queue_front_iter
:
    public std::iterator<std::output_iterator_tag, void, void, void, void>
{
  public:
    typedef typename Queue::value_type value_type;

    queue_front_iter(Queue& q) : q_(&q) { }
    queue_front_iter() : q_(static_cast<Queue*>(NULL)) { }

    queue_front_iter& operator *() { return *this; }
    queue_front_iter& operator ++() { return *this; }
    queue_front_iter& operator ++(int) { return *this; }
    queue_front_iter& operator =(const value_type& value);

    bool operator ==(const queue_front_iter& y) { return q_ == y.q_; }
    bool operator !=(const queue_front_iter& y) { return q_ != y.q_; }

  private:
    Queue* q_;
};

template <typename Queue>
class queue_back_iter
:
    public std::iterator<std::input_iterator_tag, void, void, void, void>
{
  public:
    typedef typename Queue::value_type value_type;

    class value
    {
      public:
        value(value_type v) : v_(v) { }
        value_type operator *() const { return v_; }
      private:
        value_type v_;
    };

    queue_back_iter(Queue& q) : q_(&q) { if ( q_ ) next(); }
    queue_back_iter() : q_(static_cast<Queue*>(NULL)) { }

    const value_type& operator *() const { return v_; }
    const value_type* operator ->() const { return &v_; }
    queue_back_iter& operator ++() { next(); return *this; }
    value operator ++(int) { value t = v_; next(); return t; }

    bool operator ==(const queue_back_iter& y)
    { return q_ == y.q_; }
    bool operator !=(const queue_back_iter& y)
    { return q_ != y.q_; }

  private:
    void next();

    Queue* q_;
    value_type v_;
};

CXX0X_ENUM_CLASS queue_op_status
{
    success = 0,
    empty,
    full,
    closed
};

#if 0
template <typename Value>
class queue_common
{
  public:
    typedef Value& reference;
    typedef const Value& const_reference;
    typedef Value value_type;

    virtual void close() = 0;
    virtual bool is_closed() = 0;
    virtual bool is_empty() = 0;

    virtual const char* name() = 0;

  protected:
    virtual ~queue_common();
};

template <typename Value>
queue_common<Value>::~queue_common() CXX0X_DEFAULTED_EASY
#endif

template <typename Queue>
class generic_queue_front
{
  public:
    typedef typename Queue::value_type value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;

    typedef queue_front_iter<generic_queue_front> iterator;
    typedef const queue_front_iter<generic_queue_front> const_iterator;

    //FIX generic_queue_front() CXX0X_DEFAULTED_EASY
    generic_queue_front(Queue& queue) : queue_(&queue) { }
    generic_queue_front(Queue* queue) : queue_(queue) { }
    generic_queue_front(const generic_queue_front& other)
        CXX0X_DEFAULTED_HARD( : queue_(other.queue_) { } )
    generic_queue_front& operator =(const generic_queue_front& other)
        CXX0X_DEFAULTED_HARD( { queue_ = other.queue_; } )

    void close() { queue_->close(); }
    bool is_closed() { return queue_->is_closed(); }
    bool is_empty() { return queue_->is_empty(); }
    const char* name() { return queue_->name(); }

    iterator begin() { return iterator(*this); }
    iterator end() { return iterator(); }
    const iterator cbegin() { return const_iterator(*this); }
    const iterator cend() { return const_iterator(); }

    void push(const value_type& x)
        { queue_->push(x); }
    queue_op_status try_push(const value_type& x)
        { return queue_->try_push(x); }
    queue_op_status wait_push(const value_type& x)
        { return queue_->wait_push(x); }
#ifdef HAS_CXX0X_RVREF
    void push(value_type&& x)
        { queue_->push( std::move(x) ); }
    queue_op_status try_push(value_type&& x)
        { return queue_->try_push( std::move(x) ); }
    queue_op_status wait_push(value_type&& x)
        { return queue_->wait_push( std::move(x) ); }
#endif

  protected:
    Queue* queue_;
};

template <typename Queue>
class generic_queue_back
{
  public:
    typedef typename Queue::value_type value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;

    typedef queue_back_iter<generic_queue_back> iterator;
    typedef queue_back_iter<generic_queue_back> const_iterator;

    //FIX generic_queue_back() CXX0X_DEFAULTED_EASY
    generic_queue_back(Queue& queue) : queue_(&queue) { }
    generic_queue_back(Queue* queue) : queue_(queue) { }
    generic_queue_back(const generic_queue_back& other)
        CXX0X_DEFAULTED_HARD( : queue_(other.queue_) { } )
    generic_queue_back& operator =(const generic_queue_back& other)
        CXX0X_DEFAULTED_HARD( { queue_ = other.queue_; } )

    void close() { queue_->close(); }
    bool is_closed() { return queue_->is_closed(); }
    bool is_empty() { return queue_->is_empty(); }
    const char* name() { return queue_->name(); }

    iterator begin() { return iterator(*this); }
    iterator end() { return iterator(); }
    const iterator cbegin() { return const_iterator(*this); }
    const iterator cend() { return const_iterator(); }

    value_type value_pop()
        { return queue_->value_pop(); }
    queue_op_status try_pop(value_type& x)
        { return queue_->try_pop(x); }
    queue_op_status wait_pop(value_type& x)
        { return queue_->wait_pop(x); }

  protected:
    Queue* queue_;
};

template <typename Queue>
queue_front_iter<Queue>&
queue_front_iter<Queue>::operator =(const value_type& value)
{
    queue_op_status s = q_->wait_push(value);
    if ( s != CXX0X_ENUM_QUAL(queue_op_status)success ) {
        q_ = NULL;
        throw s;
    }
    return *this;
}

template <typename Queue>
void
queue_back_iter<Queue>::next()
{
    queue_op_status s = q_->wait_pop(v_);
    if ( s == CXX0X_ENUM_QUAL(queue_op_status)closed )
        q_ = NULL;
}

template <typename Value>
class queue_base
{
  public:
    typedef Value value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;

    typedef queue_back_iter<queue_base> iterator;
    typedef queue_back_iter<queue_base> const_iterator;

    virtual ~queue_base() { }

    iterator begin() { return queue_back_iter<queue_base>(*this); }
    iterator end() { return queue_back_iter<queue_base>(); }
    const iterator cbegin() { return queue_back_iter<queue_base>(*this); }
    const iterator cend() { return queue_back_iter<queue_base>(); }

    virtual void close() = 0;
    virtual bool is_closed() = 0;
    virtual bool is_empty() = 0;

    virtual const char* name() = 0;

    virtual void push(const Value& x) = 0;
    virtual queue_op_status try_push(const Value& x) = 0;
    virtual queue_op_status wait_push(const Value& x) = 0;
#ifdef HAS_CXX0X_RVREF
    virtual void push(Value&& x) = 0;
    virtual queue_op_status try_push(Value&& x) = 0;
    virtual queue_op_status wait_push(Value&& x) = 0;
#endif

    virtual Value value_pop() = 0;
    virtual queue_op_status try_pop(Value&) = 0;
    virtual queue_op_status wait_pop(Value&) = 0;
};

//TODO(crowl): Use template aliases for queue_front and queue_back?

template <typename Value>
class queue_front
: public generic_queue_front< queue_base<Value> >
{
  public:
    queue_front() CXX0X_DEFAULTED_EASY
    queue_front(queue_base<Value>& queue)
        : generic_queue_front< queue_base<Value> >(queue) { }
    queue_front(queue_base<Value>* queue)
        : generic_queue_front< queue_base<Value> >(queue) { }
    queue_front(const queue_front<Value>& other)
        : generic_queue_front< queue_base<Value> >(other.queue_) { }
};

template <typename Value>
class queue_back
: public generic_queue_back< queue_base<Value> >
{
  public:
    queue_back() CXX0X_DEFAULTED_EASY
    queue_back(queue_base<Value>& queue)
        : generic_queue_back< queue_base<Value> >(queue) { }
    queue_back(queue_base<Value>* queue)
        : generic_queue_back< queue_base<Value> >(queue) { }
    queue_back(const queue_back<Value>& other)
        : generic_queue_back< queue_base<Value> >(other.queue_) { }
};

template <typename Queue>
class queue_wrapper
:
    public virtual queue_base <typename Queue::value_type>
{
    Queue* ptr;

  public:
    typedef typename Queue::value_type value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;

    queue_wrapper(Queue* arg)
    : ptr(arg)
    { }

    queue_wrapper(Queue& arg)
    : ptr(&arg)
    { }

    virtual ~queue_wrapper()
    { }

    virtual void close()
    { ptr->close(); }

    virtual bool is_closed()
    { return ptr->is_closed(); }

    virtual bool is_empty()
    { return ptr->is_empty(); }

    virtual const char* name()
    { return ptr->name(); }

    virtual void push(const value_type& x)
    { ptr->push(x); }

    virtual queue_op_status try_push(const value_type& x)
    { return ptr->try_push(x); }

    virtual queue_op_status wait_push(const value_type& x)
    { return ptr->wait_push(x); }

#ifdef HAS_CXX0X_RVREF
    virtual void push(value_type&& x)
    { ptr->push(x); }

    virtual queue_op_status try_push(value_type&& x)
    { return ptr->try_push(x); }

    virtual queue_op_status wait_push(value_type&& x)
    { return ptr->wait_push(x); }
#endif

    virtual value_type value_pop()
    { return ptr->value_pop(); }

    virtual queue_op_status try_pop(value_type& x)
    { return ptr->try_pop(x); }

    virtual queue_op_status wait_pop(value_type& x)
    { return ptr->wait_pop(x); }

    queue_front<value_type> front()
    { return queue_front<value_type>(this); }

    queue_back<value_type> back()
    { return queue_back<value_type>(this); }
};

template <typename Value>
class queue_counted
:
    public queue_base<Value>
{
  public:
    queue_counted() : f_(0), r_(0) { }
    virtual ~queue_counted() { }

    void inc_front() { f_++; }
    void inc_back() { r_++; }
    bool dec_front() { return --f_ == 0; }
    bool dec_back() { return --r_ == 0; }
    bool no_front() { return f_ == 0; }
    bool no_back() { return r_ == 0; }

  private:
    std::atomic<int> f_;
    std::atomic<int> r_;
};

template <typename Value>
class shared_queue_front
{
  public:
    typedef Value value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;

    typedef queue_front_iter<shared_queue_front> iterator;
    typedef const queue_front_iter<shared_queue_front> const_iterator;

    //FIX shared_queue_front()
    //FIX     : queue_(NULL) { }
    shared_queue_front(queue_counted<value_type>* queue)
        : queue_(queue) { queue->inc_front(); }
    shared_queue_front(const shared_queue_front& other)
        : queue_(other.queue_) { queue_->inc_front(); }
#ifdef HAS_CXX0X_RVREF
    shared_queue_front(shared_queue_front&& other)
        : queue_(other.queue_) { other.queue_ = NULL; }
#endif

  private:
    void release()
    {
        if ( queue_ != NULL && queue_->dec_front() ) {
            queue_->close();
            if ( queue_->no_back() ) {
                delete queue_;
            }
        }
    }

  public:
    ~shared_queue_front() { release(); }

    shared_queue_front& operator =(const shared_queue_front& other)
    {
        if ( this != &other ) {
            release();
            queue_ = other->queue_;
            if ( queue_ != NULL )
                queue_->inc_front();
        }
        return *this;
    }
#ifdef HAS_CXX0X_RVREF
    shared_queue_front& operator =(shared_queue_front&& other)
    {
        if ( this != &other ) {
            release();
            queue_ = other->queue_;
            other->queue_ == NULL;
        }
        return *this;
    }
#endif

    void close() { queue_->close(); }
    bool is_closed() { return queue_->is_closed(); }
    bool is_empty() { return queue_->is_empty(); }
    const char* name() { return queue_->name(); }

    iterator begin() { return iterator(*this); }
    iterator end() { return iterator(); }
    const iterator cbegin() { return const_iterator(*this); }
    const iterator cend() { return const_iterator(); }

    void push(const value_type& x)
        { queue_->push(x); }
    queue_op_status try_push(const value_type& x)
        { return queue_->try_push(x); }
    queue_op_status wait_push(const value_type& x)
        { return queue_->wait_push(x); }
#ifdef HAS_CXX0X_RVREF
    void push(value_type&& x)
        { queue_->push( std::move(x) ); }
    queue_op_status try_push(value_type&& x)
        { return queue_->try_push( std::move(x) ); }
    queue_op_status wait_push(value_type&& x)
        { return queue_->wait_push( std::move(x) ); }
#endif

  private:
    queue_counted<value_type>* queue_;
};

template <typename Value>
class shared_queue_back
{
  public:
    typedef Value value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;

    typedef queue_back_iter<shared_queue_back> iterator;
    typedef queue_back_iter<shared_queue_back> const_iterator;

    //FIX shared_queue_back()
    //FIX     : queue_(NULL) { }
    shared_queue_back(queue_counted<value_type>* queue)
        : queue_(queue) { queue->inc_back(); }
    shared_queue_back(const shared_queue_back& other)
        : queue_(other.queue_) { queue_->inc_back(); }
#ifdef HAS_CXX0X_RVREF
    shared_queue_back(shared_queue_back&& other)
        : queue_(other.queue_) { other.queue_ = NULL; }
#endif

  private:
    void release()
    {
        if ( queue_ != NULL && queue_->dec_back() ) {
            queue_->close();
            if ( queue_->no_front() ) {
                delete queue_;
            }
        }
    }

  public:
    ~shared_queue_back() { release(); }

    shared_queue_back& operator =(const shared_queue_back& other)
    {
        if ( this != &other ) {
            release();
            queue_ = other->queue_;
            if ( queue_ != NULL )
                queue_->inc_back();
        }
        return *this;
    }
#ifdef HAS_CXX0X_RVREF
    shared_queue_back& operator =(shared_queue_back&& other)
    {
        if ( this != &other ) {
            release();
            queue_ = other->queue_;
            other->queue_ = NULL;
        }
        return *this;
    }
#endif

    void close() { queue_->close(); }
    bool is_closed() { return queue_->is_closed(); }
    bool is_empty() { return queue_->is_empty(); }
    const char* name() { return queue_->name(); }

    iterator begin() { return iterator(*this); }
    iterator end() { return iterator(); }
    const iterator cbegin() { return const_iterator(*this); }
    const iterator cend() { return const_iterator(); }

    value_type value_pop()
        { return queue_->value_pop(); }
    queue_op_status try_pop(value_type& x)
        { return queue_->try_pop(x); }
    queue_op_status wait_pop(value_type& x)
        { return queue_->wait_pop(x); }

  private:
    queue_counted<value_type>* queue_;
};

template <typename Queue>
class queue_owner
:
    public queue_counted <typename Queue::value_type>
{
    Queue* ptr;

  public:
    typedef typename Queue::value_type value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;

    queue_owner(const queue_owner&) CXX0X_DELETED
    queue_owner(Queue* arg) : ptr(arg) { }

    virtual ~queue_owner() { delete ptr; }

    queue_front<value_type> front()
        { return queue_front<value_type>(this); }
    queue_back<value_type> back()
        { return queue_back<value_type>(this); }

    virtual void close() { ptr->close(); }
    virtual bool is_closed() { return ptr->is_closed(); }
    virtual bool is_empty() { return ptr->is_empty(); }
    virtual const char* name() { return ptr->name(); }

    virtual void push(const value_type& x)
        { ptr->push(x); }
    virtual queue_op_status try_push(const value_type& x)
        { return ptr->try_push(x); }
    virtual queue_op_status wait_push(const value_type& x)
        { return ptr->wait_push(x); }

#ifdef HAS_CXX0X_RVREF
    virtual void push(value_type&& x)
        { ptr->push(x); }
    virtual queue_op_status try_push(value_type&& x)
        { return ptr->try_push(x); }
    virtual queue_op_status wait_push(value_type&& x)
        { return ptr->wait_push(x); }
#endif

    virtual value_type value_pop()
        { return ptr->value_pop(); }
    virtual queue_op_status try_pop(value_type& x)
        { return ptr->try_pop(x); }
    virtual queue_op_status wait_pop(value_type& x)
        { return ptr->wait_pop(x); }
};


template <typename Queue>
class queue_object
:
    public queue_counted <typename Queue::value_type>
{
    Queue obj_;

  public:
    typedef typename Queue::value_type value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;

    queue_object(const queue_object&) CXX0X_DELETED
#ifdef HAS_CXX0X_VARIADIC_TMPL
    template <typename ... Args>
    queue_object(Args ... args) : obj_(args...) { }
#else
    template <typename Arg>
    queue_object(Arg arg) : obj_(arg) { }
    template <typename Arg1, typename Arg2>
    queue_object(Arg1 arg1, Arg2 arg2) : obj_(arg1, arg2) { }
#endif

    virtual ~queue_object() { }

    operator queue_front<value_type>() //TODO(crowl) Really?
        { return queue_front<value_type>(this); }
    queue_front<value_type> front()
        { return queue_front<value_type>(this); }
    queue_back<value_type> back()
        { return queue_back<value_type>(this); }

    virtual void close() { obj_.close(); }
    virtual bool is_closed() { return obj_.is_closed(); }
    virtual bool is_empty() { return obj_.is_empty(); }
    virtual const char* name() { return obj_.name(); }

    virtual void push(const value_type& x)
        { obj_.push(x); }
    virtual queue_op_status try_push(const value_type& x)
        { return obj_.try_push(x); }
    virtual queue_op_status wait_push(const value_type& x)
        { return obj_.wait_push(x); }

#ifdef HAS_CXX0X_RVREF
    virtual void push(value_type&& x)
        { obj_.push(x); }
    virtual queue_op_status try_push(value_type&& x)
        { return obj_.try_push(x); }
    virtual queue_op_status wait_push(value_type&& x)
        { return obj_.wait_push(x); }
#endif

    virtual value_type value_pop()
        { return obj_.value_pop(); }
    virtual queue_op_status try_pop(value_type& x)
        { return obj_.try_pop(x); }
    virtual queue_op_status wait_pop(value_type& x)
        { return obj_.wait_pop(x); }
};


#ifdef HAS_CXX0X_VARIADIC_TMPL

template <typename Queue, typename ... Args>
std::pair< shared_queue_front<typename Queue::value_type>,
           shared_queue_back<typename Queue::value_type> >
share_queue_ends(Args ... args)
{
  typedef typename Queue::value_type elemtype;
  CXX0X_AUTO_VAR( q, new queue_object<Queue>(args...) );
  return std::make_pair(shared_queue_front<elemtype>(q),
                        shared_queue_back<elemtype>(q));
}

#else

template <typename Queue, typename Arg>
std::pair< shared_queue_front<typename Queue::value_type>,
           shared_queue_back<typename Queue::value_type> >
share_queue_ends(Arg arg)
{
  typedef typename Queue::value_type elemtype;
  CXX0X_AUTO_VAR( q, new queue_object<Queue>(arg) );
  return std::make_pair(shared_queue_front<elemtype>(q),
                        shared_queue_back<elemtype>(q));
}

template <typename Queue, typename Arg1, typename Arg2>
std::pair< shared_queue_front<typename Queue::value_type>,
           shared_queue_back<typename Queue::value_type> >
share_queue_ends(Arg1 arg1, Arg2 arg2)
{
  typedef typename Queue::value_type elemtype;
  CXX0X_AUTO_VAR( q, new queue_object<Queue>(arg1, arg2) );
  return std::make_pair(shared_queue_front<elemtype>(q),
                        shared_queue_back<elemtype>(q));
}

#endif

} // namespace gcl

#endif
