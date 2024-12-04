//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __KMER_VECTOR_HPP__
#define __KMER_VECTOR_HPP__

#include "array_vector.hpp"
#include "config.hpp"

#ifdef SPADES_USE_JEMALLOC

# include <jemalloc/jemalloc.h>

#endif

namespace adt {

template<class Seq>
class KMerVector {
private:
    typedef typename Seq::DataType ElTy;

    ElTy *realloc() {
        storage_ = (ElTy*)::realloc(storage_, sizeof(ElTy) * capacity_ * el_sz_);
        return storage_;
    }

public:
    typedef typename array_vector<ElTy>::reference reference;
    typedef typename array_vector<ElTy>::value_type value_type;
    typedef typename array_vector<ElTy>::iterator iterator;
    typedef typename array_vector<ElTy>::const_iterator const_iterator;

    typedef array_less<ElTy> less2_fast;
    typedef array_equal_to<ElTy> equal_to;

    explicit KMerVector(unsigned K, size_t capacity = 1)
            : K_(K), size_(0), capacity_(std::max(capacity, (size_t) 1)), el_sz_(Seq::GetDataSize(K)),
              storage_(NULL),
              vector_(realloc(), size_, el_sz_) {
    }

    KMerVector(KMerVector &&that)
            : K_(that.K_), size_(that.size_), capacity_(that.capacity_), el_sz_(that.el_sz_),
              storage_(that.storage_),
              vector_(storage_, size_, el_sz_) {
        that.storage_ = NULL;
    }

    KMerVector(const KMerVector &that)
            : K_(that.K_), size_(that.size_), capacity_(that.capacity_), el_sz_(that.el_sz_),
              storage_(NULL),
              vector_(realloc(), size_, el_sz_) {
        memcpy(storage_, that.storage_, size_ * sizeof(ElTy) * el_sz_);
    }

    ~KMerVector() {
        free(storage_);
    }

    KMerVector &operator=(const KMerVector &that) {
        if (this != &that) {
            K_ = that.K_;
            size_ = that.size_;
            capacity_ = that.capacity_;
            el_sz_ = that.el_sz_;

            storage_ = NULL;
            realloc();
            memcpy(storage_, that.storage_, size_ * sizeof(ElTy) * el_sz_);

            vector_.set_data(storage_);
            vector_.set_size(size_);
        }

        return *this;
    }

    void push_back(const ElTy *data) {
        if (capacity_ == size_)
            reserve(capacity_ * 2);

        vector_[size_] = data;
        size_ += 1;
        vector_.set_size(size_);
    }

    void push_back(const Seq &s) {
        push_back(s.data());
    }

    void push_back(reference s) {
        push_back(s.data());
    }

    void push_back(const value_type &s) {
        push_back(s.data());
    }

    void reserve(size_t amount) {
        if (capacity_ < amount) {
            capacity_ = amount;
            vector_.set_data(realloc());
        }
    }

    void clear() {
        size_ = 0;
        vector_.set_size(size_);
    }

    void shrink(size_t sz) {
        size_ = sz;
        vector_.set_size(size_);
    }

    void shrink_to_fit() {
        capacity_ = std::max(size_, size_t(1));
        vector_.set_data(realloc());
    }

    iterator begin() {
        return vector_.begin();
    }

    const_iterator begin() const {
        return vector_.begin();
    }

    iterator end() {
        return vector_.end();
    }

    const_iterator end() const {
        return vector_.end();
    }

    reference back() {
        return vector_[size_-1];
    }
    
    ElTy *data() {
        return storage_;
    }

    const ElTy *data() const {
        return storage_;
    }

    size_t size() const {
        return size_;
    }

    size_t el_size() const {
        return el_sz_;
    }

    size_t el_data_size() const {
        return el_sz_ * sizeof(ElTy);
    }

    size_t capacity() const {
        return capacity_;
    }

    const ElTy *operator[](size_t idx) const {
        return vector_[idx];
    }

    unsigned K() const {
        return K_;
    }

private:
    unsigned K_;
    size_t size_;
    size_t capacity_;
    size_t el_sz_;
    ElTy *storage_;
    array_vector<ElTy> vector_;
};

} //adt
#endif /* __KMER_VECTOR_HPP */
