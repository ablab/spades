//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * kmer_set.hpp
 *
 *  Created on: Jul 19, 2012
 *      Author: alex
 */

#ifndef KMER_SET_HPP_
#define KMER_SET_HPP_


#include "runtime_k.hpp"
#include "kmer_hash_vector.hpp"


namespace runtime_k  {


class IKmerSet {

public:

    typedef RtSeq input_value_type;

    virtual ~IKmerSet() {

    }

    virtual IKmerSet * copy() const = 0;

    virtual bool empty() const = 0;

    virtual size_t size() const = 0;

    virtual size_t count(const input_value_type& kmer_seq) const = 0;

    virtual bool insert(const input_value_type& val) = 0;

    virtual void transfer(IKmerHashVector * vec, size_t thread_num) = 0;

    virtual void clear() = 0;

    virtual void erase() = 0;

    virtual void to_file(const std::string& s) const = 0;

    virtual size_t get_k() const = 0;

    virtual bool contains(const input_value_type& val) = 0;

};



class KmerSet {

public:

    typedef IKmerSet base_set_type;

private:

    base_set_type * data_;

public:

    typedef KmerSet set_type;

    typedef base_set_type::input_value_type input_value_type;

    KmerSet(size_t k);

    KmerSet(base_set_type * set): data_(set) {
    }

    KmerSet(const set_type& set) {
        data_ = set.data_->copy();
    }

    set_type& operator=(const set_type& set) {
        if (set.data_ != data_) {
            delete data_;
            data_ = set.data_->copy();
        }

        return *this;
    }

    ~KmerSet() {
       delete data_;
    }


    bool empty() const {
        return data_->empty();
    }

    size_t size() const {
        return data_->size();
    }

    size_t count(const input_value_type& kmer_seq) const {
        return data_->count(kmer_seq);
    }

    bool insert(const input_value_type& val) {
        return data_->insert(val);
    }

    void transfer(const KmerHashVector& vec, size_t thread_num) {
        data_->transfer(vec.get_data(), thread_num);
    }

    void clear() {
        data_->clear();
    }

    void erase() {
        data_->erase();
    }

    void to_file(const std::string& s) const {
        data_->to_file(s);
    }

    size_t get_k() const {
        return data_->get_k();
    }

    base_set_type * get_data() const {
        return data_;
    }

    bool contains(const input_value_type& val) {
    	return data_->contains(val);
    }

};



// ================================= SET IMPLEMENTATION =================================



template <size_t size_>
class KmerSetImpl: public IKmerSet {

public:

    typedef TypeContainerImpl<size_> type_container;

    typedef typename type_container::set_type set_type;

    typedef IKmerSet base_type;

    typedef typename base_type::input_value_type input_value_type;

private:

    set_type data_;

    size_t k_;

public:

    KmerSetImpl(size_t k, size_t n): data_(n), k_(k) {
    }

    KmerSetImpl(size_t k): data_(), k_(k) {
    }

    virtual base_type * copy() const {
        return new KmerSetImpl<size_>(*this);
    }

    virtual bool empty() const {
        return data_.empty();
    }

    virtual size_t size() const {
        return data_.size();
    }

    virtual size_t count(const input_value_type& kmer_seq) const {
        return data_.count(type_container::from_sequence(kmer_seq));
    }

    virtual bool insert(const input_value_type& val) {
        return data_.insert(type_container::from_sequence(val)).second;
    }

    virtual void transfer(IKmerHashVector * vec, size_t thread_num) {
        VERIFY_MSG(vec->get_k() == k_, "Unable to transfer vector to set of different k values");

        //KmerHashVectorImpl<size_> * vec_impl = (KmerHashVectorImpl<size_> *) vec;
        KmerHashVectorImpl<size_> * vec_impl = dynamic_cast< KmerHashVectorImpl<size_> *>(vec);
        data_.insert(vec_impl->get_data()[thread_num].begin(), vec_impl->get_data()[thread_num].end());

//        for (auto iter = vec_impl->get_data()[thread_num].begin(); iter != vec_impl->get_data()[thread_num].end(); ++iter) {
//            data_.insert(*iter);
//        }
    }

    virtual void clear() {
        data_.clear();
    }

    virtual void erase() {
        data_.erase(data_.begin(), data_.end());
    }

    virtual size_t get_k() const {
        return k_;
    }

    virtual void to_file(const std::string& s) const {
        ofstream kmeros;
        kmeros.open(s.c_str());
        for (auto iter = data_.begin(); iter != data_.end(); ++iter) {
            kmeros << *iter << std::endl;
        }
        kmeros.close();
    }

    virtual bool contains(const input_value_type& val) {
    	return (data_.find(type_container::from_sequence(val)) != data_.end());
    }

    const set_type& get_data() const {
        return data_;
    }

};

// ================================= SET FACTORIES =================================
// Single factory interface
class SingleKmerSetFactory {

public:

    virtual IKmerSet * GetSet(size_t k, size_t capacity) const = 0;

    virtual IKmerSet * GetSet(size_t k) const = 0;

    virtual ~SingleKmerSetFactory() {

    }

};


// Single factory for specific k and value
template <size_t ts_>
class SingleKmerSetFactoryImpl: public SingleKmerSetFactory {

public:

    virtual IKmerSet * GetSet(size_t k, size_t capacity) const {
        VERIFY_MSG(GET_UPPER_BOUND(k) == GET_K_BY_TS(ts_), k << " -> " << GET_UPPER_BOUND(k) << ", " << ts_ << " -> " << GET_K_BY_TS(ts_));

        return new KmerSetImpl< GET_K_BY_TS(ts_) >(k, capacity);
    }

    virtual IKmerSet * GetSet(size_t k) const {
        VERIFY_MSG(GET_UPPER_BOUND(k) == GET_K_BY_TS(ts_), k << " -> " << GET_UPPER_BOUND(k) << ", " << ts_ << " -> " << GET_K_BY_TS(ts_));

        return new KmerSetImpl< GET_K_BY_TS(ts_) >(k);
    }

};

//Factory genetator
template<size_t ts_>
class SetGenerator {

public:

    static void GenerateSets(std::vector< SingleKmerSetFactory* > & factories) {
        factories[ts_] = new SingleKmerSetFactoryImpl<ts_>();
        SetGenerator<ts_ - 1> :: GenerateSets (factories);
    }
};

//Terminating factory generator
template<>
class SetGenerator<MIN_TS> {

public:

    static void GenerateSets(std::vector< SingleKmerSetFactory* > & factories) {
        factories[MIN_TS] = new SingleKmerSetFactoryImpl<MIN_TS>();
    }
};


//Lazy singleton for factory for every required value
class KmerSetFactory {

private:

    std::vector < SingleKmerSetFactory* > single_factories_;

    KmerSetFactory() {
        VERIFY_MSG(MIN_K <= MAX_K, "Invalid K value range");

        single_factories_ = std::vector < SingleKmerSetFactory* >(MAX_TS + 1);
        SetGenerator<MAX_TS>::GenerateSets(single_factories_);
    }

public:

    static KmerSetFactory& GetInstance() {
        static KmerSetFactory instance;

        return instance;
    }

    KmerSet GetSet(size_t k, size_t capacity) {
        VERIFY_MSG(k >= MIN_K && k <= MAX_K, "K value " + ToString(k) + " is not supported, should be >= " +
                ToString(MIN_K) + " and <= " + ToString(MAX_K));

        return KmerSet(single_factories_[GET_T_ELEMENTS_NUMBER(k)]->GetSet(k, capacity));
    }

    KmerSet GetSet(size_t k) {
        VERIFY_MSG(k >= MIN_K && k <= MAX_K, "K value " + ToString(k) + " is not supported, should be >= " +
                ToString(MIN_K) + " and <= " + ToString(MAX_K));

        return KmerSet(single_factories_[GET_T_ELEMENTS_NUMBER(k)]->GetSet(k));
    }

    IKmerSet * GetRawSet(size_t k) {
        VERIFY_MSG(k >= MIN_K && k <= MAX_K, "K value " + ToString(k) + " is not supported, should be >= " +
                ToString(MIN_K) + " and <= " + ToString(MAX_K));

        return single_factories_[GET_T_ELEMENTS_NUMBER(k)]->GetSet(k);
    }
};


// Main set getters
KmerSet GetSet(size_t k, size_t capacity) {
    return KmerSetFactory::GetInstance().GetSet(k, capacity);
}

KmerSet GetSet(size_t k) {
    return KmerSetFactory::GetInstance().GetSet(k);
}

KmerSet::KmerSet(size_t k): data_(KmerSetFactory::GetInstance().GetRawSet(k)) {
}

} /* namespace runtime_k */


#endif /* KMER_SET_HPP_ */
