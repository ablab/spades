#pragma once
//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "perfect_hash_map.hpp"
#include "adt/cqf.hpp"

namespace utils {
template<class K, class traits = kmers::kmer_index_traits<K>, class StoringType = SimpleStoring>
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
    void BinRead(Reader &reader) {
        KeyBase::BinRead(reader);
    }

  private:
    std::unique_ptr<ValueStorage> values_;

    friend struct CQFHashMapBuilder;
};

}
