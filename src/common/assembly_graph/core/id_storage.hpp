//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "id_distributor.hpp"
#include "utils/verify.hpp"

#include <atomic>
#include <cstdlib>

namespace omnigraph {

template<class T>
class IdStorage {
  private:
    void resize(size_t N) {
        VERIFY(N > storage_size_);
        T *new_storage = (T*)malloc(N * sizeof(T));
        // Note: we cannot iterate over id's here as resize() could
        // be called during create() and therefore we might end with
        // uninitialized newly created id as well.
        for (uint64_t id = bias_; id < storage_size_; ++id) {
            if (!id_distributor_.occupied(id))
                continue;

            T *val = &storage_[id];
            ::new(new_storage + id) T(std::move(*val));
        }

        if (storage_)
            free(storage_);
        storage_ = new_storage;
        storage_size_ = N;
    }

  public:
    typedef omnigraph::ReclaimingIdDistributor::id_iterator id_iterator;
    typedef T value_type;

    IdStorage(uint64_t bias)
            : size_(0), bias_(bias), storage_(nullptr), storage_size_(0), id_distributor_(bias) {
        resize(id_distributor_.size() + bias_);
    }

    ~IdStorage() {
        if (storage_) {
            for (uint64_t id : id_distributor_.ids()) {
                T *val = &storage_[id];
                val->~T();
            }

            free(storage_);
        }
    }

    id_iterator id_begin() const { return id_distributor_.begin(); }
    id_iterator id_end() const { return id_distributor_.end(); }
    uint64_t max_id() const { return id_distributor_.max_id(); }

    void reserve(size_t sz) {
        if (storage_size_ >= sz + bias_)
            return;

        id_distributor_.resize(sz);
        resize(sz + bias_);
    }

    size_t size() const noexcept { return size_; }

    bool contains(uint64_t id) const {
        return id < storage_size_ && id_distributor_.occupied(id);
    }

    template<typename... ArgTypes>
    uint64_t create(ArgTypes &&... args) {
        uint64_t id = id_distributor_.allocate();

        while (storage_size_ < id + 1)
            resize(storage_size_ * 2 + 1);

        new(storage_ + id) T(std::forward<ArgTypes>(args)...);;
        size_ += 1;

        return id;
    }

    template<typename... ArgTypes>
    uint64_t emplace(uint64_t at, ArgTypes &&... args) {
        // One MUST call reserve before using emplace()
        VERIFY(!id_distributor_.occupied(at));

        id_distributor_.acquire(at);
        new(storage_ + at) T(std::forward<ArgTypes>(args)...);;
        size_.fetch_add(1);

        return at;
    }

    void erase(uint64_t id) {
        T *v = &storage_[id];

        v->~T();

        id_distributor_.release(id);
        size_ -= 1;
    }

    T& at(uint64_t id) const noexcept {
        return storage_[id];
    }

    uint64_t reserved() const { return id_distributor_.size(); }
    void clear_state() { id_distributor_.clear_state(); }

  private:
    std::atomic<size_t> size_;
    uint64_t bias_;
    T *storage_;
    size_t storage_size_;
    omnigraph::ReclaimingIdDistributor id_distributor_;
};

}
