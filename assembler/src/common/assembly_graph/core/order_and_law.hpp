//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <boost/utility.hpp>

#include <ostream>
#include <unordered_set>
#include <unordered_map>
#include "utils/stacktrace.hpp"
#include <algorithm>
#include <map>
#include "utils/openmp_wrapper.h"
#include "folly/PackedSyncPtr.h"


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

    IdSegmentStorage() : min_value_(0), size_(0) { }

private:
    IdSegmentStorage(size_t min_value, size_t size) : min_value_(min_value), size_(size) { }

    size_t min_value_;
    size_t size_;
};

// Id distributor for pure_pointer. Singleton.
class LocalIdDistributor : public IdDistributor, boost::noncopyable {
    friend class PeriodicIdDistributor;

    static const size_t INITIAL_MAX_INT_ID = 2;
public:
    size_t GetId() {
        return max_int_id_++;
    }

    IdSegmentStorage Reserve(size_t size) {
        max_int_id_ += size;
        return IdSegmentStorage(max_int_id_ - size, size);
    }

    IdSegmentStorage ReserveUpTo(size_t max) {
        VERIFY(max_int_id_ == INITIAL_MAX_INT_ID);
        max_int_id_ = max;
        return IdSegmentStorage(0, max);
    }

//  static GlobalIdDistributor &GetInstance() {
//    static GlobalIdDistributor instance(INITIAL_MAX_INT_ID);
//    return instance;
//  }

    size_t GetMax() const {
        return max_int_id_;
    }

    LocalIdDistributor(size_t min_id_value = INITIAL_MAX_INT_ID) : max_int_id_(min_id_value) { }

private:
    size_t max_int_id_;
};

/* id distributor used for concurrent algorithms.
* each thread use their own PeriodicIdDistributor with period equals to
* the quantity of threads. After thread's job is done Synchronize call are required
* to increase id in GlobalIdDistributor.
*/
class PeriodicIdDistributor : public IdDistributor {

public:
    PeriodicIdDistributor(LocalIdDistributor &id_distributor, size_t first_id, size_t period)
            : id_distributor_(id_distributor), cur_id_(first_id), period_(period) {
    }

    virtual size_t GetId() {
        size_t id = cur_id_;
        cur_id_ += period_;

        return id;
    }

    void Synchronize() const {
        size_t &global_max_id = id_distributor_.max_int_id_;
        global_max_id = std::max(cur_id_, global_max_id);
    }

private:
    LocalIdDistributor &id_distributor_;
    size_t cur_id_;
    size_t period_;
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

    explicit pure_pointer(T *ptr)
            : int_id_(size_t(ptr)) {
        ptr_.init(ptr, MAX_THREAD_CNT);
        VERIFY(int_id_ < 2);
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

//template<class T>
//struct Comparator
//{
//  typedef pure_pointer<T> pointer_type_t;
//
//  bool operator()(pointer_type_t const& a, pointer_type_t const& b) const {
//    return a.get() < b.get();
//  }
//};

template<class T>
struct Hash {
    typedef pure_pointer<T> pointer_type_t;
    std::hash<T *> inner_hash_;

    size_t operator()(pointer_type_t const &a) const {
        return inner_hash_(a.get());
    }
};

template<class It>
struct iterator_wrapper {
    typedef typename It::value_type value_type;
    typedef typename It::difference_type difference_type;
    typedef typename It::reference reference;
    typedef typename It::pointer pointer;

    explicit iterator_wrapper(It it) : it_(it) { }

    reference   operator*() const { return it_.operator*(); }

    pointer   operator->() const { return it_.operator->(); }

    bool operator==(const iterator_wrapper &rhs) const { return it_ == rhs.it_; }

    bool operator!=(const iterator_wrapper &rhs) const { return it_ != rhs.it_; }

private:
    It it_;
};

template<class T>
struct set {
    typedef Hash<typename T::type> hash_t;
    typedef std::unordered_set<T, hash_t> base_set_t;
    typedef typename base_set_t::value_type value_type;

    typedef iterator_wrapper<typename base_set_t::iterator> iterator;
    typedef iterator_wrapper<typename base_set_t::const_iterator> const_iterator;

public:
    set() : base_set_(10, hash_t()) {
    }

    template<class It>
    set(It begin, It end) : base_set_(begin, end, 10, hash_t()) {
    }

    const_iterator begin() const { return const_iterator(base_set_.begin()); }

    const_iterator end() const { return const_iterator(base_set_.end()); }

    iterator begin() { return iterator(base_set_.begin()); }

    iterator end() { return iterator(base_set_.end()); }

    const_iterator find(const T &key) const { return const_iterator(base_set_.find(key)); }

    iterator find(const T &key) { return iterator(base_set_.find(key)); }

    size_t count(T const &item) const { return base_set_.count(item); }

    std::pair<iterator, bool> insert(value_type const &item) {
        const std::pair<iterator, bool> &ret = base_set_.insert(item);
        return make_pair(iterator(ret.first), ret.second);
    }

    template<class It>
    void insert(It first, It last) { base_set_.insert(first, last); }

    size_t erase(const T &x) { return base_set_.erase(x); }

    void clear() { base_set_.clear(); }

    size_t size() const { return base_set_.size(); }

    bool operator==(const set &rhs) const {
        if (this->size() != rhs.size())
            return false;

        for (auto i = base_set_.begin(), j = rhs.base_set_.begin();
             i != base_set_.end() && j != rhs.base_set_.end();
             ++i, ++j) {
            if (*i != *j)
                return false;
        }

        return true;
    }

    bool operator!=(const set &rhs) const {
        return !(*this == rhs);
    }

    template<class Comparator>
    void Copy(std::set<T, Comparator> &container) const {
        container.insert(base_set_.begin(), base_set_.end());
    }

private:
    base_set_t base_set_;
};


template<class Key, class Value>
struct map {
    typedef Hash<typename Key::type> hash_t;
    typedef std::unordered_map<Key, Value, hash_t> base_map_t;
    typedef typename base_map_t::value_type value_type;

    typedef iterator_wrapper<typename base_map_t::iterator> iterator;
    typedef iterator_wrapper<typename base_map_t::const_iterator> const_iterator;

public:
    map()
            : base_map_(10, hash_t()) {
    }

    template<class It>
    map(It begin, It end)
            : base_map_(begin, end, 10, hash_t()) {
    }

    const_iterator begin() const { return const_iterator(base_map_.begin()); }

    const_iterator end() const { return const_iterator(base_map_.end()); }

    iterator begin() { return iterator(base_map_.begin()); }

    iterator end() { return iterator(base_map_.end()); }

    const_iterator find(const Key &key) const {
        return const_iterator(base_map_.find(key));
    }

    iterator find(const Key &key) { return iterator(base_map_.find(key)); }

    size_t count(Key const &item) const { return base_map_.count(item); }

    Value &operator[](Key const &x) { return base_map_[x]; }

    std::pair<iterator, bool> insert(value_type const &value) {
        std::pair<iterator, bool> ret = base_map_.insert(value);
        return make_pair(iterator(ret.first), ret.second);
    }

    template<class It>
    void insert(It first, It last) { base_map_.insert(first, last); }

    size_t erase(Key const &x) { return base_map_.erase(x); }

    void clear() { base_map_.clear(); }

    size_t size() const { return base_map_.size(); }

    bool operator==(const map &rhs) const {
        if (size() != rhs.size())
            return false;

        for (auto i = base_map_.begin(), j = rhs.base_map_.begin();
             i != base_map_.end() && j != rhs.base_map_.end();
             ++i, ++j) {
            if (*i != *j)
                return false;
        }

        return true;
    }

    bool operator!=(const map &rhs) const {
        return !(*this == rhs);
    }

    template<class Comparator>
    void Copy(std::map<Key, Value, Comparator> &container) const {
        container.insert(base_map_.begin(), base_map_.end());
    }

private:
    base_map_t base_map_;
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

template<class T, class Comparator>
class PairComparator {
private:
    Comparator comparator_;
public:
    PairComparator(Comparator comparator) : comparator_(comparator) {
    }

    bool operator()(std::pair<T, T> a, std::pair<T, T> b) const {
        return a.first == b.first ? comparator_(a.second, b.second) : comparator_(a.first, b.first);
    }
};

//
//template<typename T, class Comparator>
//class MixedComparator {
//private:
//  Comparator c1_;
//  Comparator c2_;
//public:
//  MixedComparator(const Comparator &c1, const Comparator &c2) : c1_(c1), c2_(c2) {
//  }
//
//  bool operator()(const T &a, const T &b) const {
//    if(c1_.IsAFAKE(a) || c1_.IsAFAKE(b)) {
//      if(c1_.IsAFAKEMin(a))
//        return !c1_.IsAFAKEMin(b);
//      if(c1_.IsAFAKEMax(b))
//        return c1_.IsAFAKEMax(a);
//      return false;
//    }
//    if(c1_.IsValidId(a) && c1_.IsValidId(b))
//      return c1_(a, b);
//    if(c1_.IsValidId(a))
//      return true;
//    if(c1_.IsValidId(b))
//      return false;
//    if(c2_.IsValidId(a) && c2_.IsValidId(b)) {
//      return c2_(a, b);
//    }
//    VERIFY(false);
//    return false;
//  }
//
//  bool IsValidId(T element) {
//    return c1_.IsValid(element) || c2_.IsValid(element);
//  }
//};

template<class Container, class Comparator>
class ContainerComparator {
private:
    Comparator comparator_;
public:
    ContainerComparator(const Comparator &comparator) : comparator_(comparator) {
    }

    bool operator()(const Container &a, const Container &b) const {
        for (auto ita = a.begin, itb = b.begin(); ita != a.end() && itb != b.end(); ++ita, ++itb) {
            if (*ita != *itb)
                return comparator_(*ita, *itb);
        }
        if (a.size() < b.size()) {
            return true;
        }
        return false;
    }

};

