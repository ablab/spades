#pragma once
//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "openmp_wrapper.h"
#include "standard.hpp"

#include "io/io_helper.hpp"

#include "mph_index/kmer_index.hpp"
#include "adt/kmer_vector.hpp"

#include "libcxx/sort.hpp"

#include "boost/bimap.hpp"

#include "kmer_splitters.hpp"
#include "key_with_hash.hpp"
#include "values.hpp"
#include "storing_traits.hpp"

#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cstdint>

namespace debruijn_graph {

template<class K, class traits>
class IndexWrapper {
	static const size_t InvalidIdx = size_t(-1);
public:
	typedef size_t IdxType;
	typedef K KeyType;
	typedef traits traits_t;
protected:
	typedef KMerIndex<traits>        KMerIndexT;
	//these fields are protected only for reduction of storage in edge indices BinWrite
	KMerIndexT index_;
private:
	std::string workdir_;
	unsigned k_;

protected:
	size_t raw_seq_idx(const typename KMerIndexT::KMerRawReference s) const {
		return index_.raw_seq_idx(s);
	}

	bool valid(const size_t idx) const {
		return idx != InvalidIdx && idx < index_.size();
	}
public:
	IndexWrapper(size_t k, const std::string &workdir) : k_((unsigned) k) {
		//fixme string literal
		workdir_ = path::make_temp_dir(workdir, "kmeridx");
	}

	~IndexWrapper() {
		path::remove_dir(workdir_);
	}

	void clear() {
		index_.clear();
	}

	unsigned k() const { return k_; }

public:
	template<class Writer>
	void BinWrite(Writer &writer) const {
		index_.serialize(writer);
	}

	template<class Reader>
	void BinRead(Reader &reader, const std::string &) {
		clear();
		index_.deserialize(reader);
	}

	const std::string &workdir() const {
		return workdir_;
	}
};

template<class K, class V, class traits, class StoringType>
class PerfectHashMap : public ValueArray<V>, public IndexWrapper<K, traits> {
public:
    typedef size_t IdxType;
    typedef K KeyType;
    typedef ValueArray<V> ValueBase;
    typedef IndexWrapper<KeyType, traits> KeyBase;
    using KeyBase::index_;
    typedef typename KeyBase::KMerIndexT KMerIndexT;
    typedef typename StoringTraits<K, KMerIndexT, StoringType>::KeyWithHash KeyWithHash;

    KeyWithHash ConstructKWH(const KeyType &key) const {
        return KeyWithHash(key, index_);
    }

    bool valid(const KeyWithHash &kwh) const {
        return KeyBase::valid(kwh.idx());
    }

    PerfectHashMap(size_t k, const std::string &workdir) : KeyBase(k, workdir) {
    }

    ~PerfectHashMap() {
    }

    void clear() {
        KeyBase::clear();
        ValueBase::clear();
    }

    const V get_value(const KeyWithHash &kwh) const {
        return StoringType::get_value(*this, kwh);
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

    template<class Writer>
    void BinWrite(Writer &writer) const {
        ValueBase::BinWrite(writer);
        KeyBase::BinWrite(writer);
    }

    template<class Reader>
    void BinRead(Reader &reader, const std::string &tmp) {
        KeyBase::BinRead(reader, tmp);
        ValueBase::BinRead(reader, tmp);
    }
//todo think more about hierarchy
protected:
    template <class KmerCounter>
    void BuildIndex(KmerCounter& counter, size_t bucket_num, size_t thread_num, bool save_final = true) {
        KMerIndexBuilder<KMerIndexT> builder(this->workdir(),
                             (unsigned) bucket_num,
                             (unsigned) thread_num);
        size_t sz = builder.BuildIndex(index_, counter, save_final);
        ValueBase::resize(sz);
    }
};


template<class K, class V, class traits, class StoringType>
class KeyStoringMap : public PerfectHashMap<K, V, traits, StoringType> {
private:
    typedef PerfectHashMap<K, V, traits, StoringType> base;

public:
    typedef traits traits_t;
    typedef K KMer;
    typedef typename base::IdxType KMerIdx;
    typedef typename traits::FinalKMerStorage::iterator kmer_iterator;
    typedef typename traits::FinalKMerStorage::const_iterator const_kmer_iterator;
    typedef typename base::KeyWithHash KeyWithHash;
    using base::ConstructKWH;

private:
    typename traits::FinalKMerStorage *kmers_;

    void SortUniqueKMers() const {
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

public:

    KeyStoringMap(size_t k, const std::string &workdir)
            : base(k, workdir),
              kmers_(NULL) {
    }

    ~KeyStoringMap() {
        delete kmers_;
    }

    KMer true_kmer(KeyWithHash kwh) const {
        VERIFY(this->valid(kwh));

        auto it = this->kmers_->begin() + kwh.idx();
        return (typename traits_t::raw_create()(this->k(), *it));
    }

    void clear() {
        base::clear();
        delete kmers_;
        kmers_ = NULL;
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
        if(!kwh.is_minimal())
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

    template<class KmerCounter>
    void BuildIndex(KmerCounter& counter, size_t bucket_num,
                    size_t thread_num) {
        base::BuildIndex(counter, bucket_num, thread_num);
        VERIFY(!kmers_);
        kmers_ = counter.GetFinalKMers();
        VERIFY(kmers_);
        SortUniqueKMers();
    }
};

template<class K, class V, class traits, class StoringType>
class KeyIteratingMap : public PerfectHashMap<K, V, traits, StoringType> {
    typedef PerfectHashMap<K, V, traits, StoringType> base;

    std::string KMersFilename_;

public:
    typedef StoringType storing_type;
    typedef typename base::traits_t traits_t;
    typedef typename base::KeyType KMer;
    typedef typename base::IdxType KMerIdx;
    using base::ConstructKWH;

public:

    KeyIteratingMap(size_t k, const std::string &workdir)
            : base(k, workdir),
              KMersFilename_("") {
    }

    ~KeyIteratingMap() {
    }

    typedef MMappedFileRecordArrayIterator<typename KMer::DataType> kmer_iterator;

    kmer_iterator kmer_begin() const {
        return kmer_iterator(this->KMersFilename_, KMer::GetDataSize(base::k()));
    }

    template<class KmerCounter>
    void BuildIndex(KmerCounter& counter, size_t bucket_num,
                    size_t thread_num) {
        base::BuildIndex(counter, bucket_num, thread_num);
        KMersFilename_ = counter.GetFinalKMersFname();
    }
};

//Seq is here for partial specialization
template <class Seq, class Index>
class DeBruijnStreamKMerIndexBuilder {

};

template<class Index>
class DeBruijnStreamKMerIndexBuilder<runtime_k::RtSeq, Index> {
 public:
    typedef Index IndexT;

    template <class Streams>
    size_t BuildIndexFromStream(IndexT &index,
                                Streams &streams,
                                io::SingleStream* contigs_stream = 0) const {
        DeBruijnReadKMerSplitter<typename Streams::ReadT, StoringTypeFilter<typename IndexT::storing_type>>
                splitter(index.workdir(), index.k(), 0, streams, contigs_stream);
        KMerDiskCounter<runtime_k::RtSeq> counter(index.workdir(), splitter);

        index.BuildIndex(counter, 16, streams.size());
        return 0;
    }
};

//fixme makes hierarchy a bit strange
template <class Index, class Seq = typename Index::KMer>
class DeBruijnGraphKMerIndexBuilder;

template <class Index>
class DeBruijnGraphKMerIndexBuilder<Index, runtime_k::RtSeq> {
 public:
  typedef Index IndexT;

  template<class Graph>
  void BuildIndexFromGraph(IndexT &index, const Graph &g, size_t read_buffer_size = 0) const {
      DeBruijnGraphKMerSplitter<Graph, StoringTypeFilter<typename Index::storing_type>> splitter(index.workdir(), index.k(),
                                                g, read_buffer_size);
      KMerDiskCounter<runtime_k::RtSeq> counter(index.workdir(), splitter);
      index.BuildIndex(counter, 16, 1);
  }
};

}
