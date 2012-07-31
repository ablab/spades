/*
 * runtime_map.hpp
 *
 *  Created on: Jun 21, 2012
 *      Author: andrey
 */

#ifndef RUNTIME_K_HPP_
#define RUNTIME_K_HPP_

#include "sequence/sequence.hpp"
#include "sequence/seq.hpp"
#include "sequence/simple_seq.hpp"
#include "sequence/rtseq.hpp"


#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "k_range.hpp"

namespace runtime_k {

#define T_SIZE sizeof(seq_element_type)

#define GET_T_ELEMENTS_NUMBER(value) ((value - 1) / (T_SIZE << 2) + 1)

#define GET_K_BY_TS(value) (value * (T_SIZE << 2))

#define GET_UPPER_BOUND(value) GET_K_BY_TS( GET_T_ELEMENTS_NUMBER(value) )


const size_t UPPER_BOUND = GET_UPPER_BOUND(MAX_K); //((MAX_K - 1) / (sizeof(seq_element_type) << 2) + 1) * (sizeof(seq_element_type) << 2);

const size_t MAX_TS = GET_T_ELEMENTS_NUMBER(MAX_K);

const size_t MIN_TS = GET_T_ELEMENTS_NUMBER(MIN_K);


typedef RuntimeSeq<UPPER_BOUND> RtSeq;


//Basic types and sequence <---> kmer functions
template <size_t size_>
class TypeContainerImpl {
public:
    typedef SimpleSeq<size_> Kmer;

    typedef unordered_set<Kmer, typename Kmer::hash, typename Kmer::equal_to> set_type;

    typedef std::vector<Kmer> vector_type;

    static Kmer from_sequence(const RtSeq& seq) {
        return seq.get_sseq<size_>();
    }

    static RtSeq to_sequence(const Kmer& kmer, size_t k = size_) {
        return RtSeq(kmer, k);
    }
};


template <size_t size_, typename Value>
class TypeValueContainerImpl: public TypeContainerImpl<size_> {

public:
    typedef TypeContainerImpl<size_> base;

    typedef typename base::Kmer Kmer;

    typedef typename base::set_type set_type;

    typedef unordered_map<Kmer, Value, typename Kmer::hash, typename Kmer::equal_to> map_type;

};




// ================================= VECTOR INTERFACE =================================
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
        data_ = vec.data_->copy();
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




// ================================= SET INTERFACE =================================
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

    set_type& operator=(const set_type& map) {
        data_ = map.data_->copy();
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




// ================================= MAP ITERATOR INTERFACE =================================

// Iterator interface
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
        iter_ = iter.iter_->copy();
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
        iter_ = iter.iter_->copy();
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

    virtual void transfer(IKmerSet * set, const Value& val) = 0;


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
template <typename Value>
class KmerMap {

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
        data_ = map.data_->copy();
        return *this;
    }

    ~KmerMap() {
       delete data_;
    }

    void transfer(const KmerSet& set, const Value& val) {
        data_->transfer(set.get_data(), val);
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


    typedef IKmerMap<Value> base_type;

    typedef typename base_type::key_type key_type;

    typedef typename base_type::value_type value_type;


    typedef KmerMapIteratorImpl<size_, Value> iterator_impl;

    typedef typename base_type::iterator_type iterator_type;

    typedef KmerConstMapIteratorImpl<size_, Value> const_iterator_impl;

    typedef typename base_type::const_iterator_type const_iterator_type;

private:

    map_type* data_;

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

    /*virtual*/ ~KmerMapImpl() {
    	delete data_;
    }

    virtual void transfer(IKmerSet * set, const Value& val) {
        VERIFY_MSG(set->get_k() == k_, "Unable to transfer set to map of different k values");

        //KmerSetImpl<size_> * set_impl = (KmerSetImpl<size_> *) set;
        KmerSetImpl<size_> * set_impl = dynamic_cast< KmerSetImpl<size_> *> (set);
        for (auto iter = set_impl->get_data().begin(); iter != set_impl->get_data().end(); ++iter) {
            data_->insert(make_pair(*iter, val));
        }
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

#endif /* RUNTIME_K_HPP_ */
