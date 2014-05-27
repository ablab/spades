//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * kmer_map.hpp
 *
 *  Created on: Jul 19, 2012
 *      Author: Alexander Opeykin
 */

#ifndef KMER_MAP_HPP_
#define KMER_MAP_HPP_


#include "runtime_k.hpp"


namespace runtime_k {

template <typename Value>
class IKmerMapIterator {

public:
    virtual ~IKmerMapIterator() {

    }

    typedef pair<const RtSeq, const Value&> value_type;

    typedef IKmerMapIterator<Value> iterator_type;


    virtual iterator_type * operator++() = 0;

    virtual iterator_type * operator++(int) = 0;

    virtual value_type operator*() = 0;

    virtual RtSeq first() = 0;

    virtual Value& second() = 0;

    virtual bool operator==(iterator_type * iter) const = 0;

    virtual bool operator!=(iterator_type * iter) const = 0;

    virtual iterator_type * copy() const = 0;

    virtual size_t get_k() const = 0;

};

//Iterator pointer wrapper
template <typename Value>
class KmerMapIterator {

public:

    typedef IKmerMapIterator<Value> base_iterator_type;

private:

    base_iterator_type * iter_;

public:

    typedef typename base_iterator_type::value_type value_type;

    typedef KmerMapIterator<Value> iterator_type;


    KmerMapIterator(base_iterator_type * iter): iter_(iter) {

    }

    KmerMapIterator(const iterator_type& iter) {
        iter_ = iter.iter_->copy();
    }

    iterator_type& operator=(const iterator_type& iter) {
      if (iter.iter_ != iter_) {
        delete iter_;
        iter_ = iter.iter_->copy();
      }

      return *this;
    }

    ~KmerMapIterator() {
       delete iter_;
    }


    iterator_type operator++() {
        return iterator_type(iter_->operator ++());
    }

    iterator_type operator++(int) {
        return iterator_type(iter_->operator ++(0));
    }

    value_type operator*() {
        return iter_->operator *();
    }

    const RtSeq first() {
        return iter_->first();
    }

    Value& second() {
        return iter_->second();
    }

    bool operator==(const iterator_type& iter) const {
        return iter_->operator ==(iter.iter_);
    }

    bool operator!=(const iterator_type& iter) const {
        return iter_->operator !=(iter.iter_);
    }


    size_t get_k() const {
        return iter_->get_k();
    }

    base_iterator_type * get_data() const {
        return iter_;
    }
};


// ================================= MAP CONST ITERATOR INTERFACE =================================

template <typename Value>
class IKmerConstMapIterator {

public:

    typedef pair<const RtSeq, const Value&> value_type;

    typedef IKmerConstMapIterator<Value> iterator_type;

    virtual ~IKmerConstMapIterator() {

    }

    virtual iterator_type * operator++() = 0;

    virtual iterator_type * operator++(int) = 0;

    virtual const value_type operator*() const = 0;

    virtual RtSeq first() const = 0;

    virtual const Value& second() const = 0;

    virtual bool operator==(iterator_type * iter) const = 0;

    virtual bool operator!=(iterator_type * iter) const = 0;

    virtual iterator_type * copy() const = 0;

    virtual size_t get_k() const = 0;

};


//Const iterator pointer wrapper
template <typename Value>
class KmerConstMapIterator {

public:

    typedef IKmerConstMapIterator<Value> base_iterator_type;

private:

    base_iterator_type * iter_;

public:

    typedef typename base_iterator_type::value_type value_type;

    typedef KmerConstMapIterator<Value> iterator_type;


    KmerConstMapIterator(base_iterator_type * iter): iter_(iter) {
    }

    KmerConstMapIterator(const iterator_type& iter) {
        iter_ = iter.iter_->copy();
    }

    iterator_type& operator=(const iterator_type& iter) {
      if (iter.iter_ != iter_) {
        delete iter_;
        iter_ = iter.iter_->copy();
      }

      return *this;
    }

    ~KmerConstMapIterator() {
       delete iter_;
    }


    iterator_type operator++() {
        return iterator_type(iter_->operator ++());
    }

    iterator_type operator++(int) {
        return iterator_type(iter_->operator ++(0));
    }

    const value_type operator*() const {
        return iter_->operator *();
    }

    RtSeq first() const {
        return iter_->first();
    }

    const Value& second() const {
        return iter_->second();
    }

    bool operator==(const iterator_type& iter) const {
        return iter_->operator ==(iter.iter_);
    }

    bool operator!=(const iterator_type& iter) const {
        return iter_->operator !=(iter.iter_);
    }

    size_t get_k() const {
        return iter_->get_k();
    }

    base_iterator_type * get_data() const {
        return iter_;
    }
};


// ================================= MAP INTERFACE =================================

template <typename Value>
class IKmerMap {

public:

    typedef RtSeq key_type;

    typedef pair<const key_type, Value> value_type;

    typedef IKmerMapIterator<Value> iterator_type;

    typedef IKmerConstMapIterator<Value> const_iterator_type;

    virtual ~IKmerMap() {

    }

    virtual IKmerMap<Value> * copy() const = 0;

    virtual bool empty() const = 0;

    virtual size_t size() const = 0;

    virtual size_t max_size() const = 0;


    virtual const_iterator_type * cbegin() const = 0;

    virtual iterator_type * begin() = 0;

    virtual const_iterator_type * cend() const = 0;

    virtual iterator_type * end() = 0;

    virtual Value& operator[](const key_type& kmer_seq) = 0;


    virtual const_iterator_type * cfind(const key_type& kmer_seq) const = 0;

    virtual iterator_type * find(const key_type& kmer_seq) = 0;

    virtual size_t count(const key_type& kmer_seq) const = 0;


    virtual pair<iterator_type *, bool> insert(const value_type& val) = 0;

    virtual size_t erase(const key_type& kmer_seq) = 0;

    //virtual iterator_type * erase(const_iterator_type * iter) = 0;

    virtual iterator_type * erase(iterator_type * iter) = 0;

    virtual void clear() = 0;


    virtual size_t bucket_count() const = 0;

    virtual size_t max_bucket_count() const = 0;

    virtual size_t bucket_size(size_t n) const = 0;

    virtual size_t bucket(const RtSeq& kmer_seq) const = 0;

    virtual float load_factor() const = 0;

    virtual float max_load_factor() const = 0;

    virtual void max_load_factor(float z) = 0;

    virtual void rehash(size_t n) = 0;

    virtual size_t get_k() const = 0;
};

//Map pointer wrapper
template <typename Value, typename Seq = RtSeq>
class KmerMap {
};

template <typename Value>
class KmerMap<Value, RtSeq> {

public:

    typedef IKmerMap<Value> base_map_type;

private:

    base_map_type * data_;

public:

    typedef KmerMap<Value> map_type;

    typedef typename base_map_type::key_type key_type;

    typedef typename base_map_type::value_type value_type;

    typedef KmerMapIterator<Value> iterator;

    typedef KmerConstMapIterator<Value> const_iterator;

    KmerMap(size_t k);

    KmerMap(base_map_type * map): data_(map) {
    }

    KmerMap(const map_type& map) {
        data_ = map.data_->copy();
    }

    map_type& operator=(const map_type& map) {
        if (map.data_ != data_) {
            delete data_;
            data_ = map.data_->copy();
        }

        return *this;
    }

    ~KmerMap() {
       delete data_;
    }

    bool empty() const {
        return data_->empty();
    }

    size_t size() const {
        return data_->size();
    }

    size_t max_size() const {
        return data_->max_size();
    }

    const_iterator begin() const {
        return const_iterator(data_->cbegin());
    }

    iterator begin() {
        return iterator(data_->begin());
    }

    const_iterator end() const {
        return const_iterator(data_->cend());
    }

    iterator end() {
        return iterator(data_->end());
    }

    Value& operator[](const RtSeq& kmer_seq) {
        return data_->operator [](kmer_seq);
    }

    const_iterator find(const RtSeq& kmer_seq) const {
        return const_iterator(data_->cfind(kmer_seq));
    }

    iterator find(const RtSeq& kmer_seq) {
        return iterator(data_->find(kmer_seq));
    }

    size_t count(const RtSeq& kmer_seq) const {
        return data_->count(kmer_seq);
    }

    pair<iterator, bool> insert(const value_type& val) {
        auto res = data_->insert(val);
        return make_pair(iterator(res.first), res.second);
    }

    size_t erase(const RtSeq& kmer_seq) {
        return data_->erase(kmer_seq);
    }

//    iterator erase(const const_iterator& iter) {
//        return iterator(data_->erase(iter.get_data()));
//    }

    iterator erase(const iterator& iter) {
        return iterator(data_->erase(iter.get_data()));
    }

    void clear() {
        data_->clear();
    }

    size_t bucket_count() const {
        return data_->bucket_count();
    }

    size_t max_bucket_count() const {
        return data_->max_bucket_count();
    }

    size_t bucket_size(size_t n) const {
        return data_->bucket_size(n);
    }

    size_t bucket(const RtSeq& kmer_seq) const {
        return data_->bucket(kmer_seq);
    }

    float load_factor() const {
        return data_->load_factor();
    }

    float max_load_factor() const {
        return data_->max_load_factor();
    }

    void max_load_factor(float z) {
        data_->max_load_factor(z);
    }

    void rehash(size_t n) {
        data_->rehash(n);
    }

    size_t get_k() const {
        return data_->get_k();
    }

    base_map_type& get_data() {
    	return *data_;
    }



};


// ================================= MAP ITERATOR IMPLEMENTATION =================================
template <size_t size_, typename Value>
class KmerMapIteratorImpl: public IKmerMapIterator<Value> {

    typedef TypeValueContainerImpl<size_, Value> type_container;

    typedef typename type_container::map_type map_type;

    typedef typename map_type::iterator map_iterator;


    typedef KmerMapIteratorImpl<size_, Value> iterator_impl;

    typedef IKmerMapIterator<Value> base_type;

    typedef typename base_type::value_type value_type;

private:
    map_iterator iter_;

    size_t k_;

public:

    KmerMapIteratorImpl(size_t k, const map_iterator& iter): iter_(iter), k_(k) {
    }

    virtual base_type * operator++() {
        return new iterator_impl(k_, ++iter_);
    }

    virtual base_type * operator++(int) {
        return new iterator_impl(k_, iter_++);
    }

    virtual value_type operator*() {
        return make_pair(type_container::to_sequence(iter_->first, k_), (*iter_).second);
    }


    virtual RtSeq first() {
        return type_container::to_sequence(iter_->first, k_);
    }

    virtual Value& second() {
        return iter_->second;
    }

    virtual bool operator==(base_type * iter) const {
        iterator_impl * it = dynamic_cast< iterator_impl * > (iter);
        return iter_ == it->iter_;
    }

    virtual bool operator!=(base_type * iter) const {
        return !operator ==(iter);
    }

    virtual base_type * copy() const {
        return new iterator_impl(k_, iter_);
    }

    virtual size_t get_k() const {
        return k_;
    }

    const map_iterator& get_data() const {
        return iter_;
    }
};


// ================================= MAP CONST ITERATOR IMPLEMENTATION =================================
template <size_t size_, typename Value>
class KmerConstMapIteratorImpl: public IKmerConstMapIterator<Value> {

    typedef TypeValueContainerImpl<size_, Value> type_container;

    typedef typename type_container::map_type map_type;

    typedef typename map_type::const_iterator map_iterator;


    typedef KmerConstMapIteratorImpl<size_, Value> iterator_impl;

    typedef IKmerConstMapIterator<Value> base_type;

    typedef typename base_type::value_type value_type;


private:
    map_iterator iter_;

    size_t k_;

public:

    KmerConstMapIteratorImpl(size_t k, const map_iterator& iter): iter_(iter), k_(k) {
    }

    virtual base_type * operator++() {
        return new iterator_impl(k_, ++iter_);
    }

    virtual base_type * operator++(int) {
        return new iterator_impl(k_, iter_++);
    }

    virtual const value_type operator*() const {
        return make_pair(type_container::to_sequence(iter_->first, k_), iter_->second);
    }


    virtual RtSeq first() const {
        return type_container::to_sequence(iter_->first, k_);
    }

    virtual const Value& second() const {
        return iter_->second;
    }

    virtual bool operator==(base_type * iter) const {
        iterator_impl * it = dynamic_cast< iterator_impl * > (iter);
        return iter_ == it->iter_;
    }

    virtual bool operator!=(base_type * iter) const {
        return !operator ==(iter);
    }

    virtual base_type * copy() const {
        return new iterator_impl(k_, iter_);
    }

    virtual size_t get_k() const {
        return k_;
    }

    const map_iterator& get_data() const {
        return iter_;
    }
};


// ================================= MAP IMPLEMENTATION =================================
template <size_t size_, typename Value>
class KmerMapImpl: public IKmerMap<Value> {

public:

    typedef TypeValueContainerImpl<size_, Value> type_container;

    typedef typename type_container::map_type map_type;

    typedef typename type_container::Kmer Kmer;

    typedef IKmerMap<Value> base_type;

    typedef typename base_type::key_type key_type;

    typedef typename base_type::value_type value_type;


    typedef KmerMapIteratorImpl<size_, Value> iterator_impl;

    typedef typename base_type::iterator_type iterator_type;

    typedef KmerConstMapIteratorImpl<size_, Value> const_iterator_impl;

    typedef typename base_type::const_iterator_type const_iterator_type;

private:

    map_type * data_;

    size_t k_;


public:

    KmerMapImpl(size_t k, size_t n) {
    	data_ = new map_type(n);
    	k_ = k;
    }

    KmerMapImpl(size_t k) {
    	data_ = new map_type();
    	k_ = k;
    }

    KmerMapImpl(const KmerMapImpl& map) {
        data_ = new map_type(*(map.data_));
        k_ = map.k_;
    }

    virtual ~KmerMapImpl() {
    	delete data_;
    }

    virtual base_type * copy() const {
        return new KmerMapImpl<size_, Value>(*this);
    }

    virtual bool empty() const {
        return data_->empty();
    }

    virtual size_t size() const {
        return data_->size();
    }

    virtual size_t max_size() const {
        return data_->max_size();
    }

    virtual const_iterator_type * cbegin() const {
        return new const_iterator_impl(k_, data_->begin());
    }

    virtual iterator_type * begin()  {
        return new iterator_impl(k_, data_->begin());
    }

    virtual const_iterator_type * cend() const {
        return new const_iterator_impl(k_, data_->end());
    }

    virtual iterator_type * end() {
        return new iterator_impl(k_, data_->end());
    }


    virtual Value& operator[](const key_type& kmer_seq) {
        return (*data_)[type_container::from_sequence(kmer_seq)];
    }

    Value& operator[](const Kmer& kmer) {
        return (*data_)[kmer];
    }


    virtual const_iterator_type * cfind(const key_type& kmer_seq) const {
        return new const_iterator_impl(k_, data_->find(type_container::from_sequence(kmer_seq)));
    }

    virtual iterator_type * find(const key_type& kmer_seq) {
        return new iterator_impl(k_, data_->find(type_container::from_sequence(kmer_seq)));
    }

    virtual size_t count(const key_type& kmer_seq) const {
        return data_->count(type_container::from_sequence(kmer_seq));
    }


    virtual pair<iterator_type *, bool> insert(const value_type& val) {
        auto res = data_->insert(make_pair(type_container::from_sequence(val.first), val.second));
        return make_pair(new iterator_impl(k_, res.first), res.second);
    }

    virtual size_t erase(const key_type& kmer_seq) {
        return data_->erase(type_container::from_sequence(kmer_seq));
    }

    virtual iterator_type * erase(iterator_type * iter) {
        VERIFY_MSG(iter->get_k() == k_, "Unable to erase by iterator of different k value");

        //iterator_impl * it = (iterator_impl *) iter;
        iterator_impl * it = dynamic_cast< iterator_impl * > (iter);
        return new iterator_impl(k_, data_->erase(it->get_data()));
    }

//    virtual iterator_type * erase(const_iterator_type * iter) {
//        VERIFY_MSG(iter->get_k() == size_, "Unable to erase by iterator of different k value");
//
//        //const_iterator_impl * it = (const_iterator_impl *) iter;
//        const_iterator_impl * it = dynamic_cast< const_iterator_impl * > (iter);
//        return new iterator_impl(data_->erase(it->get_data()));
//    }


    virtual void clear() {
    	delete data_;
    	data_ = new map_type();
    }

    virtual size_t bucket_count() const {
        return data_->bucket_count();
    }

    virtual size_t max_bucket_count() const {
        return data_->max_bucket_count();
    }

    virtual size_t bucket_size(size_t n) const {
        return data_->bucket_size(n);
    }

    virtual size_t bucket(const RtSeq& kmer_seq) const {
        return data_->bucket(type_container::from_sequence(kmer_seq));
    }

    virtual float load_factor() const {
        return data_->load_factor();
    }

    virtual float max_load_factor() const {
        return data_->max_load_factor();
    }

    virtual void max_load_factor(float z) {
        data_->max_load_factor(z);
    }

    virtual void rehash(size_t n) {
        data_->rehash(n);
    }

    virtual size_t get_k() const {
        return k_;
    }

};


// ================================= MAP FACTORIES =================================

// Single factory interface
template<class Value>
class SingleKmerMapFactory {

public:

    virtual IKmerMap<Value> * GetMap(size_t k, size_t capacity) const = 0;

    virtual IKmerMap<Value> * GetMap(size_t k) const = 0;

    virtual ~SingleKmerMapFactory() {

    }

};


// Single factory for specific k and value
template <size_t ts_, class Value>
class SingleKmerMapFactoryImpl: public SingleKmerMapFactory<Value> {

public:

    virtual IKmerMap<Value> * GetMap(size_t k, size_t capacity) const {
        VERIFY_MSG(GET_UPPER_BOUND(k) == GET_K_BY_TS(ts_), k << " -> " << GET_UPPER_BOUND(k) << ", " << ts_ << " -> " << GET_K_BY_TS(ts_));

        return new KmerMapImpl<GET_K_BY_TS(ts_), Value>(k, capacity);
    }

    virtual IKmerMap<Value> * GetMap(size_t k) const {
        VERIFY_MSG(GET_UPPER_BOUND(k) == GET_K_BY_TS(ts_), k << " -> " << GET_UPPER_BOUND(k) << ", " << ts_ << " -> " << GET_K_BY_TS(ts_));

        return new KmerMapImpl<GET_K_BY_TS(ts_), Value>(k);
    }

};

//Factory genetator
template<size_t ts_, class Value>
class MapGenerator {

public:

    static void GenerateMaps(std::vector< SingleKmerMapFactory<Value>* > & factories) {
        factories.at(ts_) = new SingleKmerMapFactoryImpl<ts_, Value>();
        MapGenerator<ts_ - 1, Value> :: GenerateMaps (factories);
    }
};

//Terminating factory generator
template<class Value>
class MapGenerator<MIN_TS, Value> {

public:

    static void GenerateMaps(std::vector< SingleKmerMapFactory<Value>* > & factories) {
        factories.at(MIN_TS) = new SingleKmerMapFactoryImpl<MIN_TS, Value>;
    }
};


//Lazy singleton for factory for every required value
template<class Value>
class KmerValueMapFactory {

private:

    std::vector < SingleKmerMapFactory<Value>* > single_factories_;

    KmerValueMapFactory() {
        VERIFY_MSG(MIN_K <= MAX_K, "Invalid K value range");

        single_factories_ = std::vector < SingleKmerMapFactory<Value>* >(MAX_TS + 1);
        MapGenerator<MAX_TS, Value>::GenerateMaps(single_factories_);
    }

  ~KmerValueMapFactory() {
    for (auto I = single_factories_.begin(), E = single_factories_.end(); I != E; ++ I)
      delete *I;
  }

public:

    static KmerValueMapFactory& GetInstance() {
        static KmerValueMapFactory instance;
        return instance;
    }

    KmerMap<Value> GetMap(size_t k, size_t capacity) {
        VERIFY_MSG(k >= MIN_K && k <= MAX_K, "K value " + ToString(k) + " is not supported, should be >= " +
                ToString(MIN_K) + " and <= " + ToString(MAX_K));

        return KmerMap<Value>(single_factories_[GET_T_ELEMENTS_NUMBER(k)]->GetMap(k, capacity));
    }

    KmerMap<Value> GetMap(size_t k) {
        VERIFY_MSG(k >= MIN_K && k <= MAX_K, "K value " + ToString(k) + " is not supported, should be >= " +
                ToString(MIN_K) + " and <= " + ToString(MAX_K));

        return KmerMap<Value>(single_factories_[GET_T_ELEMENTS_NUMBER(k)]->GetMap(k));
    }

    IKmerMap<Value> * GetRawMap(size_t k) {
        VERIFY_MSG(k >= MIN_K && k <= MAX_K, "K value " + ToString(k) + " is not supported, should be >= " +
                ToString(MIN_K) + " and <= " + ToString(MAX_K));

        return single_factories_[GET_T_ELEMENTS_NUMBER(k)]->GetMap(k);
    }
};


// Main map getter
template<class Value>
KmerMap<Value> GetMap(size_t k, size_t capacity) {
    return KmerValueMapFactory<Value>::GetInstance().GetMap(k, capacity);
}

template<class Value>
KmerMap<Value> GetMap(size_t k) {
    return KmerValueMapFactory<Value>::GetInstance().GetMap(k);
}

template<class Value>
KmerMap<Value>::KmerMap(size_t k): data_(KmerValueMapFactory<Value>::GetInstance().GetRawMap(k)) {
}


} /* namespace runtime_k */


#endif /* KMER_MAP_HPP_ */
