//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/stacktrace.hpp"
#include "utils/verify.hpp"
#include "utils/parallel/openmp_wrapper.h"
#include "folly/PackedSyncPtr.h"

#include <boost/utility.hpp>

#include <ostream>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <map>

namespace restricted {

//todo discuss with Anton
static const uint16_t MAX_THREAD_CNT = 128;

class IdDistributor {
public:
    virtual size_t GetId() = 0;

    virtual ~IdDistributor() {
    }
};

template<class Iter>
class ListIdDistributor : public IdDistributor {
    friend class IdSegmentStorage;

private:
    Iter left_;
    Iter right_;
    size_t shift_;
    size_t max_;

    ListIdDistributor(Iter left, Iter right, size_t shift = 0, size_t max = size_t(-1)) : left_(left),
                                                                                          right_(right),
                                                                                          shift_(shift), max_(max) {
    }

public:
    bool valid() {
        return left_ < right_;
    }

    size_t GetId() {
        size_t result = *(left_);
        VERIFY(result < max_);
        ++left_;
        return shift_ + result;
    }
};

class SegmentIterator {
private:
    size_t value_;
public:
    SegmentIterator(size_t value) : value_(value) {
    }

    size_t operator*() const {
        return value_;
    }

    void operator++() {
        value_++;
    }

    void operator++(int) {
        ++value_;
    }

    bool operator==(const SegmentIterator &that) const {
        return value_ == that.value_;
    }

    bool operator!=(const SegmentIterator &that) const {
        return value_ != that.value_;
    }
};

class IdSegmentStorage {
    friend class LocalIdDistributor;

public:
    ListIdDistributor<SegmentIterator> GetSegmentIdDistributor(size_t left, size_t right) {
        VERIFY(left < right);
        VERIFY(right <= size_);
        return ListIdDistributor<SegmentIterator>(SegmentIterator(left), SegmentIterator(right), min_value_, size_);
    }

    template<class Iter>
    ListIdDistributor<Iter> GetSegmentIdDistributor(Iter left, Iter right) {
        VERIFY(left < right);
        return ListIdDistributor<Iter>(left, right, min_value_, size_);
    }

private:
    IdSegmentStorage(size_t min_value, size_t size) : min_value_(min_value), size_(size) { }

    const size_t min_value_;
    const size_t size_;
};

// Id distributor for pure_pointer. Singleton.
class LocalIdDistributor : public IdDistributor, boost::noncopyable {
public:
    size_t GetId() {
        return next_int_id_++;
    }

    //force_zero_shift should only be used for loading from saves
    IdSegmentStorage Reserve(size_t size, bool force_zero_shift = false) {
        if (force_zero_shift)
            next_int_id_ = 0;
        next_int_id_ += size;
        return IdSegmentStorage(next_int_id_ - size, size);
    }

    size_t GetMax() const {
        return next_int_id_;
    }

    LocalIdDistributor() : next_int_id_(/*min int id*/2) { }

private:
    size_t next_int_id_;
};

template<class PurePtrT>
class PurePtrLock;

template<class PurePtrT>
class PurePtrMarker;

//todo maybe make it extend folly::PackedSyncPtr<T>?
template<class T>
struct pure_pointer {
    typedef T type;
    typedef T *pointer_type;

    explicit pure_pointer()
            : int_id_(0) {
        ptr_.init(pointer_type(0), MAX_THREAD_CNT);
    }

    explicit pure_pointer(T *ptr, IdDistributor &idDistributor)
            : int_id_(generate_id(ptr, idDistributor)) {
        ptr_.init(ptr, MAX_THREAD_CNT);
    }

//    lock_pointer_type& get_lockable() {
//        return ptr_;
//    }

    T *get() const {
        return ptr_.get();
    }

    T &operator*() const {
        return *ptr_;
    }

    T *operator->() const {
        return ptr_.get();
    }

    bool operator==(const pure_pointer &rhs) const {
        if (int_id_ == rhs.int_id_) {
            VERIFY(ptr_.get() == rhs.ptr_.get());
            return true;
        }
        return false;
    }

    bool operator!=(const pure_pointer &rhs) const {
        return !operator==(rhs);
    }

    bool operator<(const pure_pointer &rhs) const {
        return this->int_id_ < rhs.int_id_;
    }

    bool operator<=(const pure_pointer &rhs) const {
        return *this < rhs || *this == rhs;
    }

    size_t hash() const {
        return this->int_id_;
    }

    size_t int_id() const {
        return int_id_;
    }

private:
    friend class PurePtrLock<pure_pointer<T>>;

    friend class PurePtrMarker<pure_pointer<T>>;

    typedef folly::PackedSyncPtr<T> lock_pointer_type;

    static size_t generate_id(T *ptr, IdDistributor &idDistributor) {
        if (ptr == 0 || ptr == (T *) 1 || ptr == (T *) (-1)) {
            return size_t(ptr);
        }

        return idDistributor.GetId();
    }

    lock_pointer_type ptr_;

    size_t int_id_;
};

template<class LockT>
class ReEnteringLock {
    LockT &lock_;
    bool reentered_;

    uint16_t locking_thread() const {
        //don't need barrier here (as folly documentation says)
        return lock_.extra();
    }

    uint16_t current_thread() const {
        return uint16_t(omp_get_thread_num());
    }

    void Lock() {
        lock_.lock();
        lock_.setExtra(current_thread());
    }

    void Unlock() {
        lock_.setExtra(MAX_THREAD_CNT);
        lock_.unlock();
    }

public:
    ReEnteringLock(LockT &lock) :
            lock_(lock),
            reentered_(false) {
        if (locking_thread() == current_thread()) {
            reentered_ = true;
        } else {
            Lock();
        }
    }

    ~ReEnteringLock() {
        if (!reentered_) {
            Unlock();
        }
    }
};

/**
* Lock that uses a pure ptr as a target.
* Be careful NOT to pass a COPY of pure ptr you want to use as locked object!
*/
template<class PurePtrT>
class PurePtrLock {
    ReEnteringLock<typename PurePtrT::lock_pointer_type> inner_lock_;

public:
    PurePtrLock(PurePtrT &pure_ptr) :
            inner_lock_(pure_ptr.ptr_) {
    }

};

/**
* Way to "mark" pure pointer without using additional memory.
* Marking/unmarking operations are atomic
* Be careful NOT to pass a COPY of pure ptr you want to mark!
* Do not use with PurePtrLocks, they use the same space for storing data...
*/
template<class PurePtrT>
class PurePtrMarker {
    typedef typename PurePtrT::lock_pointer_type LockWithData;

    void ChangeMark(PurePtrT &pure_ptr, uint16_t new_mark) const {
        LockWithData &lock_with_data = pure_ptr.ptr_;
        lock_with_data.lock();
        lock_with_data.setExtra(new_mark);
        lock_with_data.unlock();
    }

public:

    void mark(PurePtrT &pure_ptr) const {
        ChangeMark(pure_ptr, 0);
    }

    void unmark(PurePtrT &pure_ptr) const {
        ChangeMark(pure_ptr, MAX_THREAD_CNT);
    }

    bool is_marked(const PurePtrT &pure_ptr) const {
        uint16_t curr_mark = pure_ptr.ptr_.extra();
        VERIFY(curr_mark == 0 || curr_mark == MAX_THREAD_CNT);
        return curr_mark == 0;
    }

};

template<class T>
std::ostream &operator<<(std::ostream &stream, const pure_pointer<T> &pointer) {
    stream << pointer.int_id();
    return stream;
}

} // namespace restricted

namespace std {

template<class T>
struct hash<restricted::pure_pointer<T>> {
    size_t operator()(const restricted::pure_pointer<T> &pointer) const {
        return pointer.hash();
    }
};

}

