#pragma once
//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/parallel/openmp_wrapper.h"

#include "key_with_hash.hpp"
#include "storing_traits.hpp"

#include "io/binary/binary.hpp"

#include "utils/kmer_mph/kmer_index.hpp"
#include "utils/verify.hpp"

#include <vector>
#include <cstdlib>
#include <cstdint>

namespace utils {

template<class K, class traits>
class IndexWrapper {
    static const size_t InvalidIdx = size_t(-1);
public:
    typedef size_t IdxType;
    typedef K KeyType;
    typedef traits traits_t;
protected:
    typedef kmers::KMerIndex<traits> KMerIndexT;
    //these fields are protected only for reduction of storage in edge indices BinWrite
    std::shared_ptr<KMerIndexT> index_ptr_;
private:
    unsigned k_;

protected:
    size_t raw_seq_idx(const typename KMerIndexT::KMerRawReference s) const {
        return index_ptr_->raw_seq_idx(s);
    }

    bool valid(const size_t idx) const {
        return idx != InvalidIdx && idx < index_ptr_->size();
    }
public:
    IndexWrapper(unsigned k)
            : index_ptr_(std::make_shared<KMerIndexT>()),
              k_(k) {}

    IndexWrapper(unsigned k, std::shared_ptr<KMerIndexT> index_ptr)
            : IndexWrapper(k) {
        index_ptr_ = index_ptr;
    }

    ~IndexWrapper() {}

    void clear() {
        index_ptr_->clear();
    }

    unsigned k() const { return k_; }

public:
    template<class Writer>
    void BinWrite(Writer &writer) const {
        index_ptr_->serialize(writer);
    }

    template<class Reader>
    void BinRead(Reader &reader) {
        clear();
        index_ptr_->deserialize(reader);
    }
};

template<class K, class V,
         class traits = kmers::kmer_index_traits<K>, class StoringType = SimpleStoring,
         class Container = std::vector<V>>
class PerfectHashMap : public IndexWrapper<K, traits> {
public:
    typedef size_t IdxType;
    typedef K KeyType;
    typedef IndexWrapper<KeyType, traits> KeyBase;
    using KeyBase::index_ptr_;
    typedef typename KeyBase::KMerIndexT KMerIndexT;
    typedef typename StoringTraits<K, KMerIndexT, StoringType>::KeyWithHash KeyWithHash;

    PerfectHashMap(unsigned k)
            : KeyBase(k) {}

    PerfectHashMap(unsigned k, std::shared_ptr<KMerIndexT> index_ptr)
            : KeyBase(k, index_ptr) {
        data_.resize(index_ptr_->size());
    }

    KeyWithHash ConstructKWH(const KeyType &key) const {
        return KeyWithHash(key, *index_ptr_);
    }

    bool valid(const KeyWithHash &kwh) const {
        return KeyBase::valid(kwh.idx());
    }

    const V get_value(const KeyWithHash &kwh) const {
        return StoringType::get_value(data_, kwh);
    }

    template<typename F>
    const V get_value(const KeyWithHash &kwh, const F& inverter) const {
        return StoringType::get_value(data_, kwh, inverter);
    }

    //Think twice or ask AntonB if you want to use it!
    V &get_raw_value_reference(const KeyWithHash &kwh) {
        return data_[kwh.idx()];
    }

    const V &get_raw_value_reference(const KeyWithHash &kwh) const {
        return data_[kwh.idx()];
    }

    void put_value(const KeyWithHash &kwh, const V &value) {
        StoringType::set_value(data_, kwh, value);
    }

    template<typename F>
    void put_value(const KeyWithHash &kwh, const V &value, const F& inverter) {
        StoringType::set_value(data_, kwh, value, inverter);
    }

    void clear() {
        KeyBase::clear();
        data_.clear();
    }

    template<class Writer>
    void BinWrite(Writer &writer) const {
        io::binary::BinWrite(writer, data_);
        KeyBase::BinWrite(writer);
    }

    template<class Reader>
    void BinRead(Reader &reader) {
        io::binary::BinRead(reader, data_);
        KeyBase::BinRead(reader);
    }

    size_t size() const {
        return data_.size();
    }

    auto value_begin() const { return data_.begin(); }
    auto value_cbegin() const { return data_.cbegin(); }
    auto value_end() const { return data_.end(); }
    auto value_cend() const { return data_.cend(); }

    friend struct PerfectHashMapBuilder;

  private:
    void resize(size_t sz) {
        data_.resize(sz);
    }

    Container data_;
};

}
