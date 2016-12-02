//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * kmer_hash_vector.hpp
 *
 *  Created on: Jul 19, 2012
 *      Author: alex
 */

#ifndef KMER_HASH_VECTOR_HPP_
#define KMER_HASH_VECTOR_HPP_


#include "sequence/runtime_k.hpp"
#include "kmer_map.hpp"


namespace runtime_k {

class IKmerHashVector {

protected:
    static const size_t LOAD_OVERHEAD = 1000;

    size_t      nthreads_;

    size_t      cell_size_;

public:
    typedef RtSeq input_value_type;

    IKmerHashVector(size_t nthreads)
        : nthreads_     (nthreads)
        , cell_size_    (LOAD_OVERHEAD) {
    }

    virtual ~IKmerHashVector() {

    }

    virtual IKmerHashVector * copy() const = 0;

    virtual void clear() = 0;

    virtual void clear(size_t i) = 0;

    virtual bool is_full() const = 0;

    virtual bool is_presisely_full() const = 0;

    virtual size_t capacity(size_t i) const = 0;

    virtual size_t size(size_t i) const = 0;


    virtual void insert(const input_value_type& value) = 0;

    virtual void reserve(size_t cell_size) = 0;


    virtual size_t get_k() const = 0;

    size_t get_threads_num() const
    {
        return nthreads_;
    }

    virtual void dump (KmerMap<int>& destination, size_t bucketNum) = 0;
};



class KmerHashVector {

public:

    typedef IKmerHashVector base_vector_type;

private:

    base_vector_type * data_;

public:

    typedef KmerHashVector vector_type;

    typedef base_vector_type::input_value_type input_value_type;


    KmerHashVector(size_t k, size_t nthreads);

    KmerHashVector(base_vector_type * vec): data_(vec) {
    }

    KmerHashVector(const vector_type& vec) {
        data_ = vec.data_->copy();
    }

    vector_type& operator=(const vector_type& vec) {
        if (vec.data_ != data_) {
            delete data_;
            data_ = vec.data_->copy();
        }

        return *this;
    }

    ~KmerHashVector() {
       delete data_;
    }



    bool is_full() const {
        return data_->is_full();
    }

    bool is_presisely_full() const {
        return data_->is_presisely_full();
    }

    size_t get_threads_num() const
    {
        return data_->get_threads_num();
    }


    void insert(const input_value_type& value) {
        data_->insert(value);
    }

    void clear() {
        data_->clear();
    }


    void clear(size_t i) {
        data_->clear(i);
    }

    size_t get_k() const {
        return data_->get_k();
    }

    size_t capacity(size_t i) const {
        return data_->capacity(i);
    }

    void reserve(size_t cell_size) {
        data_->reserve(cell_size);
    }

    base_vector_type * get_data() const {
        return data_;
    }

    void print_sizes() {
        for (size_t i = 0; i < data_->get_threads_num(); ++i) {
            INFO("Size " << i << ": " << data_->size(i));
        }
    }

    void dump (KmerMap<int>& destination, size_t bucketNum) {
        data_->dump(destination, bucketNum);
    }
};


// ================================= VECTOR IMPLEMENTATION =================================

template <size_t size_>
class KmerHashVectorImpl: public IKmerHashVector {

public:

    typedef TypeContainerImpl<size_> type_container;

    typedef typename type_container::Kmer Kmer;

    typedef typename type_container::vector_type vector_type;

    typedef std::vector<vector_type> data_type;

    typedef IKmerHashVector base_type;

    typedef typename base_type::input_value_type input_value_type;

private:

    data_type data_;

    size_t k_;

public:

    KmerHashVectorImpl(size_t k, size_t nthreads):
        IKmerHashVector(nthreads)
        , data_      (nthreads)
        , k_         (k)   {
    }

    virtual base_type * copy() const {
        return new KmerHashVectorImpl<size_>(*this);
    }

    virtual bool is_full() const {
        return data_[0].size() >= cell_size_;
    }

    virtual bool is_presisely_full() const {
        for (size_t i = 0; i < nthreads_; ++i) {
            if (data_[i].size() >= cell_size_)
                return true;
        }
        return false;
    }

    virtual void insert(const input_value_type& value) {
        Kmer kmer = type_container::from_sequence(value);
        data_[kmer.GetHash() % nthreads_].push_back(kmer);
    }

    virtual void clear() {
        for (size_t i = 0; i < nthreads_; ++i) {
            data_[i].clear();
        }
    }

    virtual void clear(size_t i) {
        data_[i].clear();
    }

    virtual size_t get_k() const {
        return k_;
    }

    virtual size_t capacity(size_t i) const {
        return data_[i].capacity();
    }

    virtual size_t size(size_t i) const {
        return data_[i].size();
    }

    virtual void reserve(size_t cell_size) {
        cell_size_ = cell_size;
        for (size_t i = 0; i < nthreads_; ++i) {
            data_[i].reserve(cell_size_ + LOAD_OVERHEAD);
        }
    }

    const data_type& get_data() const {
        return data_;
    }

    virtual void dump (KmerMap<int>& destination, size_t bucketNum) {
        KmerMapImpl<size_, int>& destImpl = dynamic_cast<KmerMapImpl<size_, int>&>(destination.get_data());

        for (auto it = data_[bucketNum].begin(), end = data_[bucketNum].end(); it != end; ++it) {
            ++destImpl[*it];
        }
    }
};


// ================================= VECTOR FACTORIES =================================
// Single factory interface
class SingleKmerHashVectorFactory {

public:

    virtual IKmerHashVector * GetHashVector(size_t k, size_t nthreads) const = 0;

    virtual ~SingleKmerHashVectorFactory() {

    }
};


// Single factory for specific k and value
template <size_t ts_>
class SingleKmerHashVectorFactoryImpl: public SingleKmerHashVectorFactory {

public:

    virtual IKmerHashVector * GetHashVector(size_t k, size_t nthreads) const {
        VERIFY_MSG(GET_UPPER_BOUND(k) == GET_K_BY_TS(ts_), k << " -> " << GET_UPPER_BOUND(k) << ", " << ts_ << " -> " << GET_K_BY_TS(ts_));
        //INFO(k << " -> " << GET_UPPER_BOUND(k) << ", " << ts_ << " -> " << GET_K_BY_TS(ts_));

        return new KmerHashVectorImpl< GET_K_BY_TS(ts_) >(k, nthreads);
    }

};

//Factory genetator
template<size_t ts_>
class HashVectorGenerator {

public:

    static void GenerateHashVectors(std::vector< SingleKmerHashVectorFactory* > & factories) {
        factories[ts_] = new SingleKmerHashVectorFactoryImpl<ts_>();
        HashVectorGenerator<ts_ - 1> :: GenerateHashVectors (factories);
    }
};

//Terminating factory generator
template<>
class HashVectorGenerator<MIN_TS> {

public:

    static void GenerateHashVectors(std::vector< SingleKmerHashVectorFactory* > & factories) {
        factories[MIN_TS] = new SingleKmerHashVectorFactoryImpl<MIN_TS>;
    }
};


//Lazy singleton for factory for every required value
class KmerHashVectorFactory {

private:

    std::vector < SingleKmerHashVectorFactory* > single_factories_;

    KmerHashVectorFactory() {
        VERIFY_MSG(MIN_K <= MAX_K, "Invalid K value range");

        single_factories_ = std::vector < SingleKmerHashVectorFactory* >(MAX_TS + 1);
        HashVectorGenerator<MAX_TS>::GenerateHashVectors(single_factories_);
    }

public:

    static KmerHashVectorFactory& GetInstance() {
        static KmerHashVectorFactory instance;

        return instance;
    }

    KmerHashVector GetHashVector(size_t k, size_t nthreads) {
        VERIFY_MSG(k >= MIN_K && k <= MAX_K, "K value " + ToString(k) + " is not supported, should be >= " +
                ToString(MIN_K) + " and <= " + ToString(MAX_K));

        return KmerHashVector(single_factories_[GET_T_ELEMENTS_NUMBER(k)]->GetHashVector(k, nthreads));
    }

    IKmerHashVector * GetRawHashVector(size_t k, size_t nthreads) {
        VERIFY_MSG(k >= MIN_K && k <= MAX_K, "K value " + ToString(k) + " is not supported, should be >= " +
                ToString(MIN_K) + " and <= " + ToString(MAX_K));

        return single_factories_[GET_T_ELEMENTS_NUMBER(k)]->GetHashVector(k, nthreads);
    }
};

KmerHashVector GetHashVector(size_t k, size_t nthreads) {
    return KmerHashVectorFactory::GetInstance().GetHashVector(k, nthreads);
}

KmerHashVector::KmerHashVector(size_t k, size_t nthreads): data_(KmerHashVectorFactory::GetInstance().GetRawHashVector(k, nthreads)) {
}

} //namespace runtime_k

#endif /* KMER_HASH_VECTOR_HPP_ */
