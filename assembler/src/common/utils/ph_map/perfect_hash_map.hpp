#pragma once
//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/parallel/openmp_wrapper.h"

#include "key_with_hash.hpp"
#include "values.hpp"
#include "storing_traits.hpp"

#include "adt/cqf.hpp"
#include "io/kmers/kmer_iterator.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "utils/kmer_mph/kmer_index.hpp"
#include "utils/logger/logger.hpp"
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
    typedef KMerIndex<traits> KMerIndexT;
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
    void BinRead(Reader &reader, const std::string &) {
        clear();
        index_ptr_->deserialize(reader);
    }
};

template<class K, class V, class traits = kmer_index_traits<K>, class StoringType = SimpleStoring>
class PerfectHashMap : public ValueArray<V>, public IndexWrapper<K, traits> {
public:
    typedef size_t IdxType;
    typedef K KeyType;
    typedef ValueArray<V> ValueBase;
    typedef IndexWrapper<KeyType, traits> KeyBase;
    using KeyBase::index_ptr_;
    typedef typename KeyBase::KMerIndexT KMerIndexT;
    typedef typename StoringTraits<K, KMerIndexT, StoringType>::KeyWithHash KeyWithHash;

    KeyWithHash ConstructKWH(const KeyType &key) const {
        return KeyWithHash(key, *index_ptr_);
    }

    bool valid(const KeyWithHash &kwh) const {
        return KeyBase::valid(kwh.idx());
    }

    PerfectHashMap(unsigned k) : KeyBase(k) {}

    PerfectHashMap(unsigned k, std::shared_ptr<KMerIndexT> index_ptr)
        : KeyBase(k, index_ptr) {
        ValueBase::resize(index_ptr_->size());
    }

    void clear() {
        KeyBase::clear();
        ValueBase::clear();
    }

    const V get_value(const KeyWithHash &kwh) const {
        return StoringType::get_value(*this, kwh);
    }

    template<typename F>
    const V get_value(const KeyWithHash &kwh, const F& inverter) const {
        return StoringType::get_value(*this, kwh, inverter);
    }

    //Think twice or ask AntonB if you want to use it!
    V &get_raw_value_reference(const KeyWithHash &kwh) {
        return ValueBase::operator[](kwh.idx());
    }

    const V &get_raw_value_reference(const KeyWithHash &kwh) const {
        return ValueBase::operator[](kwh.idx());
    }

    void put_value(const KeyWithHash &kwh, const V &value) {
        StoringType::set_value(*this, kwh, value);
    }

    template<typename F>
    void put_value(const KeyWithHash &kwh, const V &value, const F& inverter) {
        StoringType::set_value(*this, kwh, value, inverter);
    }

    template<class Writer>
    void BinWrite(Writer &writer) const {
        KeyBase::BinWrite(writer);
        ValueBase::BinWrite(writer);
    }

    template<class Reader>
    void BinRead(Reader &reader, const std::string &tmp) {
        KeyBase::BinRead(reader, tmp);
        ValueBase::BinRead(reader, tmp);
    }

    friend struct PerfectHashMapBuilder;
};

template<class K, class traits = kmer_index_traits<K>, class StoringType = SimpleStoring>
class CQFHashMap : public IndexWrapper<K, traits> {
public:
    typedef size_t IdxType;
    typedef K KeyType;
    typedef IndexWrapper<KeyType, traits> KeyBase;
    using KeyBase::index_ptr_;
    typedef typename KeyBase::KMerIndexT KMerIndexT;
    typedef typename StoringTraits<K, KMerIndexT, StoringType>::KeyWithHash KeyWithHash;
    typedef qf::cqf ValueStorage;

    KeyWithHash ConstructKWH(const KeyType &key) const {
        return KeyWithHash(key, *index_ptr_);
    }

    bool valid(const KeyWithHash &kwh) const {
        return KeyBase::valid(kwh.idx());
    }

    CQFHashMap(unsigned k)
            : KeyBase(k) {}

    CQFHashMap(unsigned k, std::shared_ptr<KMerIndexT> index_ptr)
            : KeyBase(k, index_ptr) {}

    ~CQFHashMap() = default;

    void clear() {
        KeyBase::clear();
    }

    uint64_t get_value(const KeyWithHash &kwh) const {
        return values_->lookup(kwh);
    }

    void add_value(const KeyWithHash &kwh, uint64_t value) {
        values_->add(kwh, value);
    }

    template<class Writer>
    void BinWrite(Writer &writer) const {
        KeyBase::BinWrite(writer);
    }

    template<class Reader>
    void BinRead(Reader &reader, const std::string &tmp) {
        KeyBase::BinRead(reader, tmp);
    }

  private:
    std::unique_ptr<ValueStorage> values_;

    friend struct CQFHashMapBuilder;
};

template<class K, class V, class traits = kmer_index_traits<K>, class StoringType = SimpleStoring>
class KeyStoringMap : public PerfectHashMap<K, V, traits, StoringType> {
private:
    typedef PerfectHashMap<K, V, traits, StoringType> base;

public:
    typedef traits traits_t;
    typedef K KMer;
    typedef typename base::IdxType KMerIdx;

    typedef MMappedRecordArrayReader<typename KMer::DataType> KMerStorage;
    typedef typename KMerStorage::iterator kmer_iterator;
    typedef typename KMerStorage::const_iterator const_kmer_iterator;

    typedef typename base::KeyWithHash KeyWithHash;
    using base::ConstructKWH;

private:
    typename traits::ResultFile kmers_file_;
    mutable std::unique_ptr<KMerStorage> kmers_;

    void SortUniqueKMers() const {
        if (!kmers_)
            kmers_.reset(new KMerStorage(*kmers_file_, KMer::GetDataSize(base::k())));

        size_t swaps = 0;
        INFO("Arranging kmers in hash map order");
        for (auto I = kmers_->begin(), E = kmers_->end(); I != E; ++I) {
            size_t cidx = I - kmers_->begin();
            size_t kidx = this->raw_seq_idx(*I);
            while (cidx != kidx) {
                auto J = kmers_->begin() + kidx;
                using std::swap;
                swap(*I, *J);
                swaps += 1;
                kidx = this->raw_seq_idx(*I);
            }
        }
        INFO("Done. Total swaps: " << swaps);
    }

protected:
    template<class Writer>
    void BinWriteKmers(Writer &writer) const {
        traits::raw_serialize(writer, this->kmers_);
    }

    template<class Reader>
    void BinReadKmers(Reader &reader, const std::string &FileName) {
        this->kmers_ = traits_t::raw_deserialize(reader, FileName);
    }

public:
    template<class Writer>
    void BinWrite(Writer &writer) const {
        base::BinWrite(writer);
        BinWriteKmers(writer);
    }

    template<class Reader>
    void BinRead(Reader &reader, const std::string &FileName) {
        base::BinRead(reader, FileName);
        BinReadKmers(reader, FileName);
    }

    KeyStoringMap(unsigned k)
            : base(k), kmers_(nullptr) {}

    KeyStoringMap(KeyStoringMap&& other)
            : base(std::move(other)), kmers_(std::move(other.kmers_)) {}

    KMer true_kmer(KeyWithHash kwh) const {
        VERIFY(this->valid(kwh));

        auto it = this->kmers_->begin() + kwh.idx();
        return (typename traits_t::raw_create()(this->k(), *it));
    }

    void clear() {
        base::clear();
        kmers_.reset(nullptr);
    }

    kmer_iterator kmer_begin() {
        return kmers_->begin();
    }
    const_kmer_iterator kmer_begin() const {
        return kmers_->cbegin();
    }

    kmer_iterator kmer_end() {
        return kmers_->end();
    }
    const_kmer_iterator kmer_end() const {
        return kmers_->cend();
    }

    bool valid(const KeyWithHash &kwh) const {
        if (!base::valid(kwh))
            return false;

        auto it = this->kmers_->begin() + kwh.idx();
        if (!kwh.is_minimal())
            return (typename traits_t::raw_equal_to()(!kwh.key(), *it));
        else
            return (typename traits_t::raw_equal_to()(kwh.key(), *it));
    }

    /**
    * Number of edges going out of the param edge's end
    */
    unsigned NextEdgeCount(const KeyWithHash &kwh) const {
        unsigned res = 0;
        for (char c = 0; c < 4; ++c)
          if (valid(kwh << c))
            res += 1;

        return res;
    }

    KeyWithHash NextEdge(const KeyWithHash &kwh) const { // returns any next edge
        for (char c = 0; c < 4; ++c) {
          if (valid(kwh << c))
            //hack for this code to work with long seqs! (oterwise return s is totally fine)
            return ConstructKWH(true_kmer(kwh));//s;
        }

        VERIFY_MSG(false, "Couldn't find requested edge!");
        return ConstructKWH(KMer(this->k()));
        // no next edges (we should request one here).
    }

    /**
    * Number of edges coming into param edge's end
    */
    unsigned RivalEdgeCount(const KeyWithHash &kwh) const {
        KeyWithHash next = kwh << 'A';
        unsigned res = 0;
        for (char c = 0; c < 4; ++c)
          if (valid(next >> c))
            res += 1;

        return res;
    }

    friend struct KeyStoringIndexBuilder;
};

template<class K, class V, class traits = kmer_index_traits<K>, class StoringType = SimpleStoring>
class KeyIteratingMap : public PerfectHashMap<K, V, traits, StoringType> {
    typedef PerfectHashMap<K, V, traits, StoringType> base;

    typename traits::ResultFile kmers_;

public:
    typedef StoringType storing_type;
    typedef typename base::traits_t traits_t;
    typedef typename base::KeyType KMer;
    typedef typename base::IdxType KMerIdx;
    using base::ConstructKWH;

public:

    KeyIteratingMap(unsigned k)
            : base(k) {}

    ~KeyIteratingMap() {}

    typedef MMappedFileRecordArrayIterator<typename KMer::DataType> kmer_iterator;

    kmer_iterator kmer_begin() const {
        VERIFY(kmers_ && "Index should be built");
        return kmer_iterator(*this->kmers_, KMer::GetDataSize(base::k()));
    }

    std::vector<kmer_iterator> kmer_begin(size_t parts) const {
        VERIFY(kmers_ && "Index should be built");
        return io::make_kmer_iterator<KMer>(*this->kmers_, base::k(), parts);
    }

    friend struct KeyIteratingIndexBuilder;
};

}
