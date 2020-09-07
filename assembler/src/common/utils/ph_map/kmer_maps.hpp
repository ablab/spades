#pragma once
//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "perfect_hash_map.hpp"
#include "io/kmers/kmer_iterator.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "utils/logger/logger.hpp"

namespace utils {

template<class K, class V, class traits = kmers::kmer_index_traits<K>, class StoringType = SimpleStoring>
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
        base::BinRead(reader);
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

template<class K, class V, class traits = kmers::kmer_index_traits<K>, class StoringType = SimpleStoring>
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
        return raw_kmer_iterator(*this->kmers_, KMer::GetDataSize(base::k()));
    }

    std::vector<kmer_iterator> kmer_begin(size_t parts) const {
        VERIFY(kmers_ && "Index should be built");
        return io::make_raw_kmer_iterator<KMer>(*this->kmers_, base::k(), parts);
    }

    friend struct KeyIteratingIndexBuilder;
};

}
